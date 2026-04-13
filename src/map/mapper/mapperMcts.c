/**CFile****************************************************************

  FileName    [mapperMcts.c]

  PackageName [MVSIS 1.3: Multi-valued logic synthesis system.]

  Synopsis    [Local MCTS refinement for technology mapping.]

  Author      [OpenAI]

***********************************************************************/

#include "mapperInt.h"
#include <stdarg.h>
#ifdef ABC_USE_OPENMP
#include <omp.h>
#endif

ABC_NAMESPACE_IMPL_START

#define MAP_MCTS_MAX_DECISIONS  10
#define MAP_MCTS_MAX_PHASECANDS 3
#define MAP_MCTS_MAX_TEMP_COMBOS ((MAP_MCTS_MAX_PHASECANDS + 1) * (MAP_MCTS_MAX_PHASECANDS + 1))
#define MAP_MCTS_MAX_COMBOS     8
#define MAP_MCTS_TFI_HOPS       2
#define MAP_MCTS_TFI_MAXNODES   128
#define MAP_MCTS_LOG_FILE       "mcts_runtime.log"

typedef struct Map_MctsNodeState_t_ Map_MctsNodeState_t;
typedef struct Map_MctsSnapshot_t_  Map_MctsSnapshot_t;
typedef struct Map_MctsCombo_t_     Map_MctsCombo_t;
typedef struct Map_MctsDecision_t_  Map_MctsDecision_t;
typedef struct Map_MctsRegion_t_    Map_MctsRegion_t;
typedef struct Map_MctsTree_t_      Map_MctsTree_t;
typedef struct Map_MctsCtx_t_       Map_MctsCtx_t;
typedef struct Map_MctsResult_t_    Map_MctsResult_t;

struct Map_MctsNodeState_t_
{
    Map_Cut_t * pCutBest[2];
    Map_Time_t  tArrival[2];
    Map_Time_t  tRequired[2];
    int         nRefAct[3];
};

struct Map_MctsSnapshot_t_
{
    int                 nNodes;
    Map_MctsNodeState_t * pStates;
    float               RequiredGlo;
    float               AreaFinal;
};

struct Map_MctsCombo_t_
{
    Map_Cut_t *         pCutBest[2];
    float               Score;
    float               Prior;
};

struct Map_MctsDecision_t_
{
    Map_Node_t *        pNode;
    int                 nCombos;
    float               PotScore;
    Map_MctsCombo_t     Combos[MAP_MCTS_MAX_COMBOS];
};

struct Map_MctsRegion_t_
{
    Map_Node_t *        pRoot;
    Map_NodeVec_t *     vNodes;
    float               Score;
    float               SharedRefs;
    float               PairSharedRefs;
    int                 nRepeatNodes;
    int                 nSharedLeaves;
    int                 nPairTfiOverlap;
    int                 nPairLeafOverlap;
    int                 nDecisions;
    Map_MctsDecision_t  Decisions[MAP_MCTS_MAX_DECISIONS];
};

struct Map_MctsTree_t_
{
    int                 Depth;
    int                 nVisits;
    float               ValueSum;
    int                 nChildren;
    float               Priors[MAP_MCTS_MAX_COMBOS];
    Map_MctsTree_t *    pChildren[MAP_MCTS_MAX_COMBOS];
};

struct Map_MctsCtx_t_
{
    Map_Man_t *         p;
    Map_MctsRegion_t *  pRegion;
    Map_MctsSnapshot_t * pSnapshot;
    Map_Node_t *        pRoot;
    float               AreaBase;
    float               DelayBound;
    float               BestValue;
    float               BestArea;
    float               BestDelay;
    int                 iRegion;
    int                 nRegions;
    int                 iSim;
    int                 BestChoices[MAP_MCTS_MAX_DECISIONS];
};

struct Map_MctsResult_t_
{
    int                 fFound;
    float               BestValue;
    float               BestArea;
    float               BestDelay;
    int                 BestChoices[MAP_MCTS_MAX_DECISIONS];
};

static FILE * s_pMapMctsLog = NULL;

static float Map_MctsDelayUsingOther( Map_Man_t * p, Map_Cut_t * pCut, int fPhase );
static int   Map_MctsTraceOpen();
static void  Map_MctsTraceClose();
static void  Map_MctsTrace( const char * pFormat, ... );
static Map_Node_t * Map_MctsRemapNode( Map_Man_t * pNew, Map_Node_t * pNode );
static Map_Man_t * Map_MctsManDup( Map_Man_t * p );
static void  Map_MctsManFree( Map_Man_t * p );
static void  Map_MctsSnapshotSave( Map_Man_t * p, Map_MctsSnapshot_t * pSnap );
static void  Map_MctsSnapshotRestore( Map_Man_t * p, Map_MctsSnapshot_t * pSnap );
static void  Map_MctsSnapshotFree( Map_MctsSnapshot_t * pSnap );
static float Map_MctsNodeSlack( Map_Node_t * pNode );
static float Map_MctsNodeArea( Map_Node_t * pNode );
static void  Map_MctsCollectTfiNode_rec( Map_Node_t * pNode, int fPhase, int nHops, Map_Node_t ** ppNodes, int * pnNodes, int nNodesMax );
static void  Map_MctsCollectTfiCut( Map_Cut_t * pCut, int fPhase, Map_Node_t ** ppNodes, int * pnNodes, int nNodesMax );
static void  Map_MctsCollectCurrentTfi( Map_Node_t * pNode, Map_Node_t ** ppNodes, int * pnNodes, int nNodesMax );
static int   Map_MctsCollectCurrentLeaves( Map_Node_t * pNode, Map_Node_t ** ppLeaves, int nLeavesMax );
static int   Map_MctsCountNodeOverlap( Map_Node_t ** ppNodes0, int nNodes0, Map_Node_t ** ppNodes1, int nNodes1, float * pSharedRefs );
static void  Map_MctsCollectAigTfiNode_rec( Map_Node_t * pNode, int nHops, Map_Node_t ** ppNodes, int * pnNodes, int nNodesMax );
static float Map_MctsNodeHotspotScore( Map_Node_t * pNode );
static void  Map_MctsScoreRegion( Map_MctsRegion_t * pRegion );
static float Map_MctsNodeReconvScore( Map_Node_t * pNode );
static float Map_MctsPhaseCutScore( Map_Cut_t * pCut, int fPhase );
static float Map_MctsComboArea( Map_Man_t * p, Map_MctsCombo_t * pCombo );
static float Map_MctsComboDelay( Map_Man_t * p, Map_MctsCombo_t * pCombo );
static void  Map_MctsPrintChoices( Map_MctsRegion_t * pRegion, int * pChoices );
static void  Map_MctsPrintDecision( Map_Man_t * p, Map_MctsDecision_t * pDecision, int iDecision );
static void  Map_MctsPrintRegion( Map_MctsCtx_t * pCtx );
static int   Map_MctsCollectSourceFanouts( Map_Man_t * p, Map_Node_t * pSource, Map_NodeVec_t * vFanouts );
static int   Map_MctsIsReconvergentSource( Map_Man_t * p, Map_Node_t * pSource );
static int   Map_MctsCollectSourceMerges( Map_Man_t * p, Map_Node_t * pSource, unsigned char * pReach, Map_NodeVec_t * vFanouts, Map_NodeVec_t * vMerges );
static void  Map_MctsCollectSourceRegion_rec( Map_Node_t * pNode, Map_Node_t * pSource, unsigned char * pReach, Map_NodeVec_t * vNodes, int nLimit );
static int   Map_MctsCompareRegions( const void * p1, const void * p2 );
static int   Map_MctsComparePhaseCuts( const void * p1, const void * p2 );
static int   Map_MctsCompareCombos( const void * p1, const void * p2 );
static void  Map_MctsSortPhaseCuts( Map_Cut_t ** ppCuts, int nCuts, int fPhase );
static int   Map_MctsCollectPhaseCuts( Map_Node_t * pNode, int fPhase, Map_Cut_t ** ppCuts, int nCutsMax );
static int   Map_MctsBuildDecision( Map_Man_t * p, Map_Node_t * pNode, Map_MctsDecision_t * pDecision );
static float Map_MctsDecisionPotentialScore( Map_MctsRegion_t * pRegion, Map_MctsDecision_t * pDecision );
static int   Map_MctsCompareDecisionScores( const void * p1, const void * p2 );
static void  Map_MctsBuildRegion( Map_Man_t * p, Map_Node_t * pRoot, Map_MctsRegion_t * pRegion );
static int   Map_MctsRecomputeCover( Map_Man_t * p );
static int   Map_MctsApplyChoices( Map_MctsRegion_t * pRegion, int * pChoices );
static int   Map_MctsCompleteChoices( Map_MctsCtx_t * pCtx, int * pChoicesIn, int * pChoicesOut );
static float Map_MctsEvaluateChoices( Map_MctsCtx_t * pCtx, int * pChoices );
static float Map_MctsSearch( Map_MctsCtx_t * pCtx, Map_MctsTree_t * pTree, int * pChoices );
static int   Map_MctsSearchOneRegion( Map_Man_t * p, Map_Node_t * pRoot, float AreaBase, float DelayBase, int iRegion, int nRegions, Map_MctsResult_t * pRes );
static void  Map_MctsTreeFree( Map_MctsTree_t * pTree );

static float Map_MctsDelayUsingOther( Map_Man_t * p, Map_Cut_t * pCut, int fPhase )    // 当某个相位没有直接实现时，利用“另一相位的 cut + 一个 inverter”来估算该相位的延迟
{
    Map_Time_t * pTime;
    if ( pCut == NULL )
        return MAP_FLOAT_LARGE;
    pTime = &pCut->M[fPhase].tArrive;
    if ( pTime == NULL || pTime->Worst >= MAP_FLOAT_LARGE/2 )
        return MAP_FLOAT_LARGE;
    if ( fPhase )
        return MAP_MAX( pTime->Fall + p->pSuperLib->tDelayInv.Rise, pTime->Rise + p->pSuperLib->tDelayInv.Fall );
    return MAP_MAX( pTime->Fall + p->pSuperLib->tDelayInv.Rise, pTime->Rise + p->pSuperLib->tDelayInv.Fall );
}

static int Map_MctsTraceOpen()
{
    if ( s_pMapMctsLog != NULL )
        return 1;
    s_pMapMctsLog = fopen( MAP_MCTS_LOG_FILE, "w" );
    return s_pMapMctsLog != NULL;
}

static void Map_MctsTraceClose()
{
    if ( s_pMapMctsLog == NULL )
        return;
    fclose( s_pMapMctsLog );
    s_pMapMctsLog = NULL;
}

static void Map_MctsTrace( const char * pFormat, ... )
{
    va_list Args;
    if ( s_pMapMctsLog == NULL )
        return;
#ifdef ABC_USE_OPENMP
#pragma omp critical(map_mcts_log)
#endif
    {
    va_start( Args, pFormat );
    vfprintf( s_pMapMctsLog, pFormat, Args );
    va_end( Args );
    fflush( s_pMapMctsLog );
    }
}

static inline const char * Map_MctsLogStr( const char * pStr )
{
    return pStr ? pStr : "<unknown>";
}

static Map_Node_t * Map_MctsRemapNode( Map_Man_t * pNew, Map_Node_t * pNode )   // 把旧映射管理器里的一个节点 pNode，重新映射到新管理器 pNew 中对应的节点，并保留“是否取反”的信息
{
    Map_Node_t * pNodeR;
    int fComp;
    if ( pNode == NULL )
        return NULL;
    pNodeR = Map_Regular( pNode );
    fComp = Map_IsComplement( pNode );
    if ( Map_NodeIsConst(pNodeR) )
        return Map_NotCond( pNew->pConst1, fComp );
    assert( pNodeR->Num >= 0 && pNodeR->Num < pNew->vMapObjs->nSize );
    return Map_NotCond( pNew->vMapObjs->pArray[pNodeR->Num], fComp );
}

static Map_Man_t * Map_MctsManDup( Map_Man_t * p )    // 深拷贝/复制整个 Map 管理器
{
    Map_Man_t * pNew;
    Map_Node_t * pNodes, * pConst1;
    Map_Node_t * pNodeOld, * pNodeNew;
    Map_Cut_t * pCutOld, * pCutNew, * pCutPrev;
    int i, k;

    pNew = ABC_ALLOC( Map_Man_t, 1 );
    *pNew = *p;
    pNew->pBins = NULL;
    pNew->vVisited = NULL;
    pNew->mmNodes = NULL;
    pNew->mmCuts = NULL;
    pNew->fUseProfile = 0;

    pConst1 = ABC_ALLOC( Map_Node_t, 1 );
    *pConst1 = *p->pConst1;
    pConst1->p = pNew;
    pConst1->pNext = NULL;
    pConst1->pCuts = NULL;
    pNew->pConst1 = pConst1;

    pNodes = ABC_ALLOC( Map_Node_t, p->vMapObjs->nSize );
    pNew->vMapObjs = Map_NodeVecAlloc( p->vMapObjs->nSize );
    for ( i = 0; i < p->vMapObjs->nSize; i++ )
    {
        pNodeOld = p->vMapObjs->pArray[i];
        pNodeNew = pNodes + i;
        *pNodeNew = *pNodeOld;
        pNodeNew->p = pNew;
        pNodeNew->pNext = NULL;
        pNodeNew->pCuts = NULL;
        pNodeNew->pCutBest[0] = NULL;
        pNodeNew->pCutBest[1] = NULL;
        Map_NodeVecPush( pNew->vMapObjs, pNodeNew );
    }
    pNew->pInputs = ABC_ALLOC( Map_Node_t *, p->nInputs );
    for ( i = 0; i < p->nInputs; i++ )
        pNew->pInputs[i] = Map_Regular( Map_MctsRemapNode( pNew, p->pInputs[i] ) );
    pNew->pOutputs = ABC_ALLOC( Map_Node_t *, p->nOutputs );
    for ( i = 0; i < p->nOutputs; i++ )
        pNew->pOutputs[i] = Map_MctsRemapNode( pNew, p->pOutputs[i] );
    pNew->vMapBufs = Map_NodeVecAlloc( Map_NodeVecReadSize(p->vMapBufs) );
    for ( i = 0; i < Map_NodeVecReadSize(p->vMapBufs); i++ )
        Map_NodeVecPush( pNew->vMapBufs, Map_Regular(Map_MctsRemapNode( pNew, Map_NodeVecReadEntry(p->vMapBufs, i) )) );

    for ( i = 0; i < pNew->vMapObjs->nSize; i++ )   // 修正每个节点内部的引用关系
    {
        pNodeNew = pNew->vMapObjs->pArray[i];
        pNodeNew->p1 = Map_MctsRemapNode( pNew, pNodeNew->p1 );
        pNodeNew->p2 = Map_MctsRemapNode( pNew, pNodeNew->p2 );
        pNodeNew->pNextE = Map_MctsRemapNode( pNew, pNodeNew->pNextE );
        pNodeNew->pRepr = Map_MctsRemapNode( pNew, pNodeNew->pRepr );
    }

    for ( i = 0; i < p->vMapObjs->nSize; i++ )      // 复制每个节点的 cut 链表
    {
        pNodeOld = p->vMapObjs->pArray[i];
        pNodeNew = pNew->vMapObjs->pArray[i];
        pCutPrev = NULL;
        for ( pCutOld = pNodeOld->pCuts; pCutOld; pCutOld = pCutOld->pNext )
        {
            pCutNew = ABC_ALLOC( Map_Cut_t, 1 );
            *pCutNew = *pCutOld;
            pCutNew->pNext = NULL;
            pCutNew->pOne = NULL;
            pCutNew->pTwo = NULL;
            for ( k = 0; k < pCutNew->nLeaves; k++ )
                pCutNew->ppLeaves[k] = Map_Regular( Map_MctsRemapNode( pNew, pCutOld->ppLeaves[k] ) );
            if ( pCutPrev == NULL )
                pNodeNew->pCuts = pCutNew;
            else
                pCutPrev->pNext = pCutNew;
            pCutPrev = pCutNew;
            if ( pNodeOld->pCutBest[0] == pCutOld )
                pNodeNew->pCutBest[0] = pCutNew;
            if ( pNodeOld->pCutBest[1] == pCutOld )
                pNodeNew->pCutBest[1] = pCutNew;
        }
    }
    return pNew;
}

static void Map_MctsManFree( Map_Man_t * p )    // 释放 Map_MctsManDup() 创建出来的整套 manager 内存
{
    Map_Cut_t * pCut, * pCutNext;
    int i;
    if ( p == NULL )
        return;
    if ( p->vMapObjs )
    {
        for ( i = 0; i < p->vMapObjs->nSize; i++ )
            for ( pCut = p->vMapObjs->pArray[i]->pCuts; pCut; pCut = pCutNext )
                pCutNext = pCut->pNext, ABC_FREE( pCut );
        if ( p->vMapObjs->nSize > 0 )
            ABC_FREE( p->vMapObjs->pArray[0] );
    }
    ABC_FREE( p->pInputs );
    ABC_FREE( p->pOutputs );
    Map_NodeVecFree( p->vMapObjs );
    Map_NodeVecFree( p->vMapBufs );
    ABC_FREE( p->pConst1 );
    ABC_FREE( p );
}

static void Map_MctsSnapshotSave( Map_Man_t * p, Map_MctsSnapshot_t * pSnap )   // 把当前 Map_Man_t 里和“映射结果/时序/引用计数”相关的运行状态保存到一个快照 pSnap 中
{
    Map_Node_t * pNode;
    int i, k;
    pSnap->nNodes  = p->vMapObjs->nSize;
    pSnap->pStates = ABC_ALLOC( Map_MctsNodeState_t, pSnap->nNodes );
    for ( i = 0; i < pSnap->nNodes; i++ )
    {
        pNode = p->vMapObjs->pArray[i];
        pSnap->pStates[i].pCutBest[0] = pNode->pCutBest[0];
        pSnap->pStates[i].pCutBest[1] = pNode->pCutBest[1];
        pSnap->pStates[i].tArrival[0] = pNode->tArrival[0];
        pSnap->pStates[i].tArrival[1] = pNode->tArrival[1];
        pSnap->pStates[i].tRequired[0] = pNode->tRequired[0];
        pSnap->pStates[i].tRequired[1] = pNode->tRequired[1];
        for ( k = 0; k < 3; k++ )
            pSnap->pStates[i].nRefAct[k] = pNode->nRefAct[k];
    }
    pSnap->RequiredGlo = p->fRequiredGlo;
    pSnap->AreaFinal   = p->AreaFinal;
}

static void Map_MctsSnapshotRestore( Map_Man_t * p, Map_MctsSnapshot_t * pSnap )    // 把快照 pSnap 里保存的节点状态和全局状态，恢复回当前 manager p
{
    Map_Node_t * pNode;
    int i, k;
    assert( pSnap->nNodes == p->vMapObjs->nSize );
    for ( i = 0; i < pSnap->nNodes; i++ )
    {
        pNode = p->vMapObjs->pArray[i];
        pNode->pCutBest[0] = pSnap->pStates[i].pCutBest[0];
        pNode->pCutBest[1] = pSnap->pStates[i].pCutBest[1];
        pNode->tArrival[0] = pSnap->pStates[i].tArrival[0];
        pNode->tArrival[1] = pSnap->pStates[i].tArrival[1];
        pNode->tRequired[0] = pSnap->pStates[i].tRequired[0];
        pNode->tRequired[1] = pSnap->pStates[i].tRequired[1];
        for ( k = 0; k < 3; k++ )
            pNode->nRefAct[k] = pSnap->pStates[i].nRefAct[k];
    }
    p->fRequiredGlo = pSnap->RequiredGlo;
    p->AreaFinal    = pSnap->AreaFinal;
}

static void Map_MctsSnapshotFree( Map_MctsSnapshot_t * pSnap )
{
    ABC_FREE( pSnap->pStates );
    pSnap->nNodes = 0;
}

static float Map_MctsNodeSlack( Map_Node_t * pNode )    // 计算一个节点当前的 slack（时序裕量）
{
    float Slack = MAP_FLOAT_LARGE;
    if ( pNode->pCutBest[0] )
        Slack = MAP_MIN( Slack, pNode->tRequired[0].Worst - pNode->tArrival[0].Worst );
    if ( pNode->pCutBest[1] )
        Slack = MAP_MIN( Slack, pNode->tRequired[1].Worst - pNode->tArrival[1].Worst );
    if ( Slack == MAP_FLOAT_LARGE )
        return 0.0;
    return Slack;
}

static float Map_MctsNodeArea( Map_Node_t * pNode )    // 估算一个节点当前实现方案对应的面积开销
{
    float Area = 0.0;
    if ( pNode->pCutBest[0] )
        Area += Map_CutGetRootArea( pNode->pCutBest[0], 0 );
    if ( pNode->pCutBest[1] )
        Area += Map_CutGetRootArea( pNode->pCutBest[1], 1 );
    if ( pNode->pCutBest[0] == NULL || pNode->pCutBest[1] == NULL )
        Area += pNode->p->pSuperLib->AreaInv;
    return Area;
}

static void Map_MctsCollectTfiNode_rec( Map_Node_t * pNode, int fPhase, int nHops, Map_Node_t ** ppNodes, int * pnNodes, int nNodesMax )    // 沿当前 cover 在 fanin 方向收集局部 TFI
{
    Map_Cut_t * pCut;
    int i, k, fLeafPhase;
    unsigned uPhaseTot;
    if ( pNode == NULL || !Map_NodeIsAnd(pNode) || pNode->pRepr )
        return;
    for ( i = 0; i < *pnNodes; i++ )
        if ( ppNodes[i] == pNode )
            return;
    if ( *pnNodes >= nNodesMax )
        return;
    ppNodes[(*pnNodes)++] = pNode;
    if ( nHops == 0 )
        return;
    pCut = pNode->pCutBest[fPhase];
    if ( pCut == NULL )
    {
        fPhase = !fPhase;
        pCut = pNode->pCutBest[fPhase];
    }
    if ( pCut == NULL )
        return;
    uPhaseTot = pCut->M[fPhase].uPhaseBest;
    for ( k = 0; k < pCut->nLeaves; k++ )
    {
        fLeafPhase = ((uPhaseTot & (1 << k)) == 0);
        if ( pCut->ppLeaves[k]->pCutBest[fLeafPhase] == NULL )
            fLeafPhase = !fLeafPhase;
        Map_MctsCollectTfiNode_rec( pCut->ppLeaves[k], fLeafPhase, nHops - 1, ppNodes, pnNodes, nNodesMax );
    }
}

static void Map_MctsCollectTfiCut( Map_Cut_t * pCut, int fPhase, Map_Node_t ** ppNodes, int * pnNodes, int nNodesMax )    // 从一个当前 best cut 出发收集其局部 TFI
{
    int i, fLeafPhase;
    unsigned uPhaseTot;
    if ( pCut == NULL )
        return;
    uPhaseTot = pCut->M[fPhase].uPhaseBest;
    for ( i = 0; i < pCut->nLeaves; i++ )
    {
        fLeafPhase = ((uPhaseTot & (1 << i)) == 0);
        if ( pCut->ppLeaves[i]->pCutBest[fLeafPhase] == NULL )
            fLeafPhase = !fLeafPhase;
        Map_MctsCollectTfiNode_rec( pCut->ppLeaves[i], fLeafPhase, MAP_MCTS_TFI_HOPS, ppNodes, pnNodes, nNodesMax );
    }
}

static void Map_MctsCollectCurrentTfi( Map_Node_t * pNode, Map_Node_t ** ppNodes, int * pnNodes, int nNodesMax )    // 收集节点当前 cover 下两个相位 best cut 的局部 TFI 并做并集
{
    if ( pNode->pCutBest[0] )
        Map_MctsCollectTfiCut( pNode->pCutBest[0], 0, ppNodes, pnNodes, nNodesMax );
    if ( pNode->pCutBest[1] )
        Map_MctsCollectTfiCut( pNode->pCutBest[1], 1, ppNodes, pnNodes, nNodesMax );
}

static int Map_MctsCollectCurrentLeaves( Map_Node_t * pNode, Map_Node_t ** ppLeaves, int nLeavesMax )    // 收集节点当前两个相位 best cut 的叶子并去重
{
    Map_Cut_t * pCut;
    int fPhase, i, k, nLeaves;
    nLeaves = 0;
    for ( fPhase = 0; fPhase < 2; fPhase++ )
    {
        pCut = pNode->pCutBest[fPhase];
        if ( pCut == NULL )
            continue;
        for ( i = 0; i < pCut->nLeaves; i++ )
        {
            for ( k = 0; k < nLeaves; k++ )
                if ( ppLeaves[k] == pCut->ppLeaves[i] )
                    break;
            if ( k < nLeaves )
                continue;
            if ( nLeaves >= nLeavesMax )
                return nLeaves;
            ppLeaves[nLeaves++] = pCut->ppLeaves[i];
        }
    }
    return nLeaves;
}

static int Map_MctsCountNodeOverlap( Map_Node_t ** ppNodes0, int nNodes0, Map_Node_t ** ppNodes1, int nNodes1, float * pSharedRefs )    // 统计两个节点集合的交集大小和共享引用数
{
    int i, k, nOverlap;
    nOverlap = 0;
    *pSharedRefs = 0.0;
    for ( i = 0; i < nNodes0; i++ )
        for ( k = 0; k < nNodes1; k++ )
            if ( ppNodes0[i] == ppNodes1[k] )
            {
                nOverlap++;
                *pSharedRefs += (float)ppNodes0[i]->nRefs;
                break;
    }
    return nOverlap;
}

static void Map_MctsCollectAigTfiNode_rec( Map_Node_t * pNode, int nHops, Map_Node_t ** ppNodes, int * pnNodes, int nNodesMax )    // 沿 AIG 结构在 fanin 方向收集局部 TFI
{
    int i;
    pNode = Map_Regular( pNode );
    if ( pNode == NULL || !Map_NodeIsAnd(pNode) || pNode->pRepr )
        return;
    for ( i = 0; i < *pnNodes; i++ )
        if ( ppNodes[i] == pNode )
            return;
    if ( *pnNodes >= nNodesMax )
        return;
    ppNodes[(*pnNodes)++] = pNode;
    if ( nHops == 0 )
        return;
    Map_MctsCollectAigTfiNode_rec( pNode->p1, nHops - 1, ppNodes, pnNodes, nNodesMax );
    Map_MctsCollectAigTfiNode_rec( pNode->p2, nHops - 1, ppNodes, pnNodes, nNodesMax );
}

static float Map_MctsNodeHotspotScore( Map_Node_t * pNode )    // 仅依据 AIG 结构中两个 fanin cone 的局部重叠，给节点打 hotspot 分
{
    Map_Node_t * pTfi0[MAP_MCTS_TFI_MAXNODES], * pTfi1[MAP_MCTS_TFI_MAXNODES];
    float SharedRefs;
    int nTfi0, nTfi1, nOverlap;
    if ( pNode == NULL || !Map_NodeIsAnd(pNode) || pNode->pRepr )
        return 0.0;
    nTfi0 = nTfi1 = 0;
    Map_MctsCollectAigTfiNode_rec( pNode->p1, MAP_MCTS_TFI_HOPS, pTfi0, &nTfi0, MAP_MCTS_TFI_MAXNODES );
    Map_MctsCollectAigTfiNode_rec( pNode->p2, MAP_MCTS_TFI_HOPS, pTfi1, &nTfi1, MAP_MCTS_TFI_MAXNODES );
    SharedRefs = 0.0;
    nOverlap = Map_MctsCountNodeOverlap( pTfi0, nTfi0, pTfi1, nTfi1, &SharedRefs );
    return 12.0f * nOverlap + SharedRefs;
}

static void Map_MctsScoreRegion( Map_MctsRegion_t * pRegion )    // 用 region 内多节点之间的 TFI/叶子交集来刻画重汇聚强度
{
    Map_Node_t * pNode0, * pNode1;
    Map_Node_t * pTfi0[MAP_MCTS_TFI_MAXNODES], * pTfi1[MAP_MCTS_TFI_MAXNODES];
    Map_Node_t * pLeaves0[2 * MAP_MCTS_MAX_PHASECANDS + 4], * pLeaves1[2 * MAP_MCTS_MAX_PHASECANDS + 4];
    float SharedRefs, ScoreCore;
    int i, k, nNodes, nTfi0, nTfi1, nLeaves0, nLeaves1;
    nNodes = Map_NodeVecReadSize( pRegion->vNodes );
    pRegion->PairSharedRefs = 0.0;
    pRegion->nPairTfiOverlap = 0;
    pRegion->nPairLeafOverlap = 0;
    for ( i = 0; i < nNodes; i++ )
    {
        pNode0 = Map_NodeVecReadEntry( pRegion->vNodes, i );
        nTfi0 = 0;
        Map_MctsCollectCurrentTfi( pNode0, pTfi0, &nTfi0, MAP_MCTS_TFI_MAXNODES );
        nLeaves0 = Map_MctsCollectCurrentLeaves( pNode0, pLeaves0, 2 * MAP_MCTS_MAX_PHASECANDS + 4 );
        for ( k = i + 1; k < nNodes; k++ )
        {
            pNode1 = Map_NodeVecReadEntry( pRegion->vNodes, k );
            nTfi1 = 0;
            Map_MctsCollectCurrentTfi( pNode1, pTfi1, &nTfi1, MAP_MCTS_TFI_MAXNODES );
            nLeaves1 = Map_MctsCollectCurrentLeaves( pNode1, pLeaves1, 2 * MAP_MCTS_MAX_PHASECANDS + 4 );
            pRegion->nPairTfiOverlap += Map_MctsCountNodeOverlap( pTfi0, nTfi0, pTfi1, nTfi1, &SharedRefs );
            pRegion->PairSharedRefs += SharedRefs;
            pRegion->nPairLeafOverlap += Map_MctsCountNodeOverlap( pLeaves0, nLeaves0, pLeaves1, nLeaves1, &SharedRefs );
        }
    }
    ScoreCore = 18.0f * pRegion->nPairTfiOverlap + 10.0f * pRegion->nPairLeafOverlap + pRegion->PairSharedRefs;
    if ( nNodes > 1 )
        ScoreCore /= (float)(nNodes - 1);
    pRegion->Score = ScoreCore + 4.0f * pRegion->nRepeatNodes + 2.0f * pRegion->nSharedLeaves + 0.25f * pRegion->SharedRefs;
}

static float Map_MctsNodeReconvScore( Map_Node_t * pNode )    // 根据当前 best cover 中两相位局部 TFI 的交集估算重汇聚/共享潜力
{
    Map_Cut_t * pCut0, * pCut1;
    Map_Node_t * pTfi0[MAP_MCTS_TFI_MAXNODES], * pTfi1[MAP_MCTS_TFI_MAXNODES];
    float SharedRefs;
    int i, k, nOverlap, nTfi0, nTfi1;
    pCut0 = pNode->pCutBest[0];
    pCut1 = pNode->pCutBest[1];
    if ( pCut0 == NULL || pCut1 == NULL )
        return 0.0;
    nTfi0 = nTfi1 = 0;
    Map_MctsCollectTfiCut( pCut0, 0, pTfi0, &nTfi0, MAP_MCTS_TFI_MAXNODES );
    Map_MctsCollectTfiCut( pCut1, 1, pTfi1, &nTfi1, MAP_MCTS_TFI_MAXNODES );
    SharedRefs = 0.0;
    nOverlap = 0;
    for ( i = 0; i < nTfi0; i++ )
        for ( k = 0; k < nTfi1; k++ )
            if ( pTfi0[i] == pTfi1[k] )
            {
                nOverlap++;
                SharedRefs += (float)pTfi0[i]->nRefs;
                break;
            }
    return 12.0f * nOverlap + SharedRefs;
}

static float Map_MctsPhaseCutScore( Map_Cut_t * pCut, int fPhase )    // 给某个 cut 在指定相位 fPhase 下算一个简单评分，用来比较 cut 的优劣：面积代价 + 很小的时序代价
{
    return pCut->M[fPhase].AreaFlow + 0.01f * pCut->M[fPhase].tArrive.Worst;
}

static float Map_MctsComboArea( Map_Man_t * p, Map_MctsCombo_t * pCombo )   // 估算一个候选组合 pCombo 的局部面积代价
{
    float Area = 0.0;
    if ( pCombo->pCutBest[0] )
        Area += Map_CutGetRootArea( pCombo->pCutBest[0], 0 );
    if ( pCombo->pCutBest[1] )
        Area += Map_CutGetRootArea( pCombo->pCutBest[1], 1 );
    if ( pCombo->pCutBest[0] == NULL || pCombo->pCutBest[1] == NULL )
        Area += p->pSuperLib->AreaInv;
    return Area;
}

static float Map_MctsComboDelay( Map_Man_t * p, Map_MctsCombo_t * pCombo )  // 估算一个候选组合 pCombo 的延迟代价，并返回两个相位中更坏的那个延迟
{
    float Delay0, Delay1;
    Delay0 = pCombo->pCutBest[0] ? pCombo->pCutBest[0]->M[0].tArrive.Worst : Map_MctsDelayUsingOther( p, pCombo->pCutBest[1], 1 );
    Delay1 = pCombo->pCutBest[1] ? pCombo->pCutBest[1]->M[1].tArrive.Worst : Map_MctsDelayUsingOther( p, pCombo->pCutBest[0], 0 );
    return MAP_MAX( Delay0, Delay1 );
}

static void Map_MctsPrintChoices( Map_MctsRegion_t * pRegion, int * pChoices )
{
    char Buffer[128];
    int nSize;
    int i;
    nSize = snprintf( Buffer, sizeof(Buffer), " choices:" );
    for ( i = 0; i < pRegion->nDecisions; i++ )
        nSize += snprintf( Buffer + nSize, sizeof(Buffer) - nSize, " n%d->%d", pRegion->Decisions[i].pNode->Num, pChoices[i] );
    snprintf( Buffer + nSize, sizeof(Buffer) - nSize, "\n" );
    Map_MctsTrace( "%s", Buffer );
}

static void Map_MctsPrintDecision( Map_Man_t * p, Map_MctsDecision_t * pDecision, int iDecision )
{
    Map_MctsCombo_t * pCombo;
    Map_Super_t * pSuper0, * pSuper1;
    const char * pName0, * pName1;
    float Area, Delay;
    int c;

    Map_MctsTrace( "    decision %d: node=%d level=%d refs=%d slack=%8.2f combos=%d pot=%8.2f\n",
        iDecision, pDecision->pNode->Num, pDecision->pNode->Level, pDecision->pNode->nRefs,
        Map_MctsNodeSlack(pDecision->pNode), pDecision->nCombos, pDecision->PotScore );
    for ( c = 0; c < pDecision->nCombos; c++ )
    {
        pCombo = &pDecision->Combos[c];
        pSuper0 = pCombo->pCutBest[0] ? pCombo->pCutBest[0]->M[0].pSuperBest : NULL;
        pSuper1 = pCombo->pCutBest[1] ? pCombo->pCutBest[1]->M[1].pSuperBest : NULL;
        pName0 = (pSuper0 && pSuper0->pRoot) ? Mio_GateReadName(pSuper0->pRoot) : "-";
        pName1 = (pSuper1 && pSuper1->pRoot) ? Mio_GateReadName(pSuper1->pRoot) : "-";
        Area = Map_MctsComboArea( p, pCombo );
        Delay = Map_MctsComboDelay( p, pCombo );
        Map_MctsTrace( "      combo %d: area=%8.2f delay=%8.2f prior=%5.2f phase0=%s/%dL phase1=%s/%dL\n",
            c, Area, Delay, pCombo->Prior, pName0, pCombo->pCutBest[0] ? pCombo->pCutBest[0]->nLeaves : 0,
            pName1, pCombo->pCutBest[1] ? pCombo->pCutBest[1]->nLeaves : 0 );
    }
}

static void Map_MctsPrintRegion( Map_MctsCtx_t * pCtx )
{
    int i;
    Map_MctsTrace( "MCTS region %d/%d: root=%d level=%d region_nodes=%d decisions=%d score=%8.2f tfi_pairs=%d leaf_pairs=%d pair_refs=%8.2f repeat=%d shared_leaves=%d shared_refs=%8.2f base_area=%11.1f delay_bound=%8.2f\n",
        pCtx->iRegion + 1, pCtx->nRegions, pCtx->pRoot->Num, pCtx->pRoot->Level,
        Map_NodeVecReadSize(pCtx->pRegion->vNodes), pCtx->pRegion->nDecisions,
        pCtx->pRegion->Score, pCtx->pRegion->nPairTfiOverlap, pCtx->pRegion->nPairLeafOverlap, pCtx->pRegion->PairSharedRefs,
        pCtx->pRegion->nRepeatNodes, pCtx->pRegion->nSharedLeaves, pCtx->pRegion->SharedRefs,
        pCtx->AreaBase, pCtx->DelayBound );
    for ( i = 0; i < pCtx->pRegion->nDecisions; i++ )
        Map_MctsPrintDecision( pCtx->p, &pCtx->pRegion->Decisions[i], i );
}

static int Map_MctsCompareRegions( const void * p1, const void * p2 )
{
    Map_MctsRegion_t * pRegion1 = (Map_MctsRegion_t *)p1;
    Map_MctsRegion_t * pRegion2 = (Map_MctsRegion_t *)p2;
    if ( pRegion1->Score < pRegion2->Score )
        return 1;
    if ( pRegion1->Score > pRegion2->Score )
        return -1;
    return pRegion2->pRoot->Level - pRegion1->pRoot->Level;
}

static int Map_MctsCollectSourceFanouts( Map_Man_t * p, Map_Node_t * pSource, Map_NodeVec_t * vFanouts )    // 收集 source 的直接内部 fanouts
{
    Map_Node_t * pNode;
    int i;
    for ( i = 0; i < p->vMapObjs->nSize; i++ )
    {
        pNode = p->vMapObjs->pArray[i];
        if ( pNode == pSource || pNode->pRepr )
            continue;
        if ( !Map_NodeIsAnd(pNode) && !Map_NodeIsBuf(pNode) )
            continue;
        if ( Map_Regular(pNode->p1) == pSource || (Map_NodeIsAnd(pNode) && Map_Regular(pNode->p2) == pSource) )
            Map_NodeVecPushUnique( vFanouts, pNode );
    }
    return Map_NodeVecReadSize( vFanouts );
}

static int Map_MctsIsReconvergentSource( Map_Man_t * p, Map_Node_t * pSource )    // 判断一个 source 的不同 fanout 分支是否会在下游重新汇合
{
    Map_NodeVec_t * vFanouts;
    Map_Node_t * pNode;
    unsigned * pTags;
    unsigned Tag0, Tag1;
    int i, RetValue;
    vFanouts = Map_NodeVecAlloc( 8 );
    pTags = ABC_CALLOC( unsigned, p->vMapObjs->nSize );
    RetValue = 0;
    if ( Map_MctsCollectSourceFanouts( p, pSource, vFanouts ) < 2 )
        goto cleanup;
    for ( i = 0; i < Map_NodeVecReadSize(vFanouts); i++ )
        pTags[Map_NodeVecReadEntry(vFanouts, i)->Num] = i + 1;
    for ( i = 0; i < p->vMapObjs->nSize; i++ )
    {
        pNode = p->vMapObjs->pArray[i];
        if ( pNode == pSource || pNode->pRepr )
            continue;
        if ( !Map_NodeIsAnd(pNode) && !Map_NodeIsBuf(pNode) )
            continue;
        Tag0 = pTags[Map_Regular(pNode->p1)->Num];
        Tag1 = Map_NodeIsAnd(pNode) ? pTags[Map_Regular(pNode->p2)->Num] : 0;
        if ( Tag0 && Tag1 && Tag0 != Tag1 )
        {
            RetValue = 1;
            break;
        }
        if ( pTags[pNode->Num] == 0 )
            pTags[pNode->Num] = Tag0 ? Tag0 : Tag1;
    }
cleanup:
    ABC_FREE( pTags );
    Map_NodeVecFree( vFanouts );
    return RetValue;
}

static int Map_MctsCollectSourceMerges( Map_Man_t * p, Map_Node_t * pSource, unsigned char * pReach, Map_NodeVec_t * vFanouts, Map_NodeVec_t * vMerges )   // 收集 source 的 merge 点并标记 source 到 merge 的可达区域
{
    Map_Node_t * pNode;
    unsigned * pTags;
    unsigned Tag0, Tag1;
    int i, nMerges;
    pTags = ABC_CALLOC( unsigned, p->vMapObjs->nSize );
    memset( pReach, 0, sizeof(unsigned char) * p->vMapObjs->nSize );
    if ( Map_MctsCollectSourceFanouts( p, pSource, vFanouts ) < 2 )
    {
        ABC_FREE( pTags );
        return 0;
    }
    for ( i = 0; i < Map_NodeVecReadSize(vFanouts); i++ )
    {
        pNode = Map_NodeVecReadEntry( vFanouts, i );
        pTags[pNode->Num] = i + 1;
        pReach[pNode->Num] = 1;
    }
    for ( i = 0; i < p->vMapObjs->nSize; i++ )
    {
        pNode = p->vMapObjs->pArray[i];
        if ( pNode == pSource || pNode->pRepr )
            continue;
        if ( !Map_NodeIsAnd(pNode) && !Map_NodeIsBuf(pNode) )
            continue;
        Tag0 = pTags[Map_Regular(pNode->p1)->Num];
        Tag1 = Map_NodeIsAnd(pNode) ? pTags[Map_Regular(pNode->p2)->Num] : 0;
        if ( Tag0 == 0 && Tag1 == 0 )
            continue;
        pReach[pNode->Num] = 1;
        if ( Tag0 && Tag1 && Tag0 != Tag1 )
        {
            Map_NodeVecPushUnique( vMerges, pNode );
            continue;
        }
        if ( pTags[pNode->Num] == 0 )
            pTags[pNode->Num] = Tag0 ? Tag0 : Tag1;
    }
    nMerges = Map_NodeVecReadSize( vMerges );
    ABC_FREE( pTags );
    return nMerges;
}

static void Map_MctsCollectSourceRegion_rec( Map_Node_t * pNode, Map_Node_t * pSource, unsigned char * pReach, Map_NodeVec_t * vNodes, int nLimit )   // 从 merge 点往回收集 source 到 merge 之间的内部节点
{
    pNode = Map_Regular( pNode );
    if ( pNode == NULL || pNode->pRepr )
        return;
    if ( pNode != pSource && !pReach[pNode->Num] )
        return;
    if ( Map_NodeIsBuf(pNode) )
    {
        Map_MctsCollectSourceRegion_rec( pNode->p1, pSource, pReach, vNodes, nLimit );
        return;
    }
    if ( !Map_NodeIsAnd(pNode) )
        return;
    if ( Map_NodeVecReadSize(vNodes) >= nLimit )
        return;
    if ( Map_NodeVecPushUnique( vNodes, pNode ) )
        return;
    if ( pNode == pSource )
        return;
    Map_MctsCollectSourceRegion_rec( pNode->p1, pSource, pReach, vNodes, nLimit );
    Map_MctsCollectSourceRegion_rec( pNode->p2, pSource, pReach, vNodes, nLimit );
}

static int Map_MctsComparePhaseCuts( const void * p1, const void * p2 )
{
    Map_Cut_t * pCut1, * pCut2;
    pCut1 = *(Map_Cut_t **)p1;
    pCut2 = *(Map_Cut_t **)p2;
    return (int)pCut1->nLeaves - (int)pCut2->nLeaves;
}

static int Map_MctsCompareCombos( const void * p1, const void * p2 )
{
    Map_MctsCombo_t * pCombo1, * pCombo2;
    pCombo1 = (Map_MctsCombo_t *)p1;
    pCombo2 = (Map_MctsCombo_t *)p2;
    if ( pCombo1->Score < pCombo2->Score )
        return -1;
    if ( pCombo1->Score > pCombo2->Score )
        return 1;
    return 0;
}

static void Map_MctsSortPhaseCuts( Map_Cut_t ** ppCuts, int nCuts, int fPhase )    // 把一组 cut 按指定相位 fPhase 的评分从小到大排序
{
    Map_Cut_t * pTemp;
    float Score1, Score2;
    int i, k;
    for ( i = 1; i < nCuts; i++ )                                                  // 分数小的在前，分数大的在后
    {
        pTemp = ppCuts[i];
        Score1 = Map_MctsPhaseCutScore( pTemp, fPhase );
        for ( k = i; k > 0; k-- )
        {
            Score2 = Map_MctsPhaseCutScore( ppCuts[k-1], fPhase );
            if ( Score2 < Score1 - MAP_FLOAT_SMALL )
                break;
            if ( Score2 < Score1 + MAP_FLOAT_SMALL && ppCuts[k-1]->nLeaves <= pTemp->nLeaves )
                break;
            ppCuts[k] = ppCuts[k-1];
        }
        ppCuts[k] = pTemp;
    }
}

static int Map_MctsCollectPhaseCuts( Map_Node_t * pNode, int fPhase, Map_Cut_t ** ppCuts, int nCutsMax )    // 为节点 pNode 在指定相位 fPhase 下收集一批可用的候选 cut，按优先级排序并返回保留数量
{
    Map_Cut_t * pCut;
    int nCuts, i;
    nCutsMax = Abc_MinInt( nCutsMax, MAP_MCTS_MAX_COMBOS );
    nCuts = 0;
    if ( pNode->pCutBest[fPhase] )
        ppCuts[nCuts++] = pNode->pCutBest[fPhase];
    for ( pCut = pNode->pCuts ? pNode->pCuts->pNext : NULL; pCut; pCut = pCut->pNext )
    {
        if ( pCut->M[fPhase].pSuperBest == NULL )
            continue;
        for ( i = 0; i < nCuts; i++ )
            if ( ppCuts[i] == pCut )
                break;
        if ( i < nCuts )
            continue;
        ppCuts[nCuts++] = pCut;
        if ( nCuts == nCutsMax )
            break;
    }
    Map_MctsSortPhaseCuts( ppCuts, nCuts, fPhase );
    return nCuts;
}

static int Map_MctsBuildDecision( Map_Man_t * p, Map_Node_t * pNode, Map_MctsDecision_t * pDecision )     // 围绕某个节点 pNode，收集正相/反相的候选 cut，枚举它们的组合并打分
{
    Map_Cut_t * pCuts0[MAP_MCTS_MAX_PHASECANDS], * pCuts1[MAP_MCTS_MAX_PHASECANDS];
    Map_MctsCombo_t Temp[MAP_MCTS_MAX_TEMP_COMBOS];
    float Delay0, Delay1, Area;
    int nCuts0, nCuts1, i, k, c, nTemp;
    int nOpt0, nOpt1;

    memset( pDecision, 0, sizeof(Map_MctsDecision_t) );
    pDecision->pNode = pNode;
    nCuts0 = Map_MctsCollectPhaseCuts( pNode, 0, pCuts0, MAP_MCTS_MAX_PHASECANDS );
    nCuts1 = Map_MctsCollectPhaseCuts( pNode, 1, pCuts1, MAP_MCTS_MAX_PHASECANDS );
    if ( nCuts0 == 0 && nCuts1 == 0 )
        return 0;

    nOpt0 = nCuts0 + (nCuts1 > 0 ? 1 : 0);
    nOpt1 = nCuts1 + (nCuts0 > 0 ? 1 : 0);
    nTemp = 0;
    for ( i = 0; i < nOpt0; i++ )
    for ( k = 0; k < nOpt1; k++ )
    {
        Temp[nTemp].pCutBest[0] = i < nCuts0 ? pCuts0[i] : NULL;
        Temp[nTemp].pCutBest[1] = k < nCuts1 ? pCuts1[k] : NULL;
        if ( Temp[nTemp].pCutBest[0] == NULL && Temp[nTemp].pCutBest[1] == NULL )
            continue;
        Area = 0.0;
        if ( Temp[nTemp].pCutBest[0] )
        {
            Area += Map_CutGetRootArea( Temp[nTemp].pCutBest[0], 0 );
            Delay0 = Temp[nTemp].pCutBest[0]->M[0].tArrive.Worst;
        }
        else
            Delay0 = Map_MctsDelayUsingOther( p, Temp[nTemp].pCutBest[1], 1 );
        if ( Temp[nTemp].pCutBest[1] )
        {
            Area += Map_CutGetRootArea( Temp[nTemp].pCutBest[1], 1 );
            Delay1 = Temp[nTemp].pCutBest[1]->M[1].tArrive.Worst;
        }
        else
            Delay1 = Map_MctsDelayUsingOther( p, Temp[nTemp].pCutBest[0], 0 );
        if ( Temp[nTemp].pCutBest[0] == NULL || Temp[nTemp].pCutBest[1] == NULL )
            Area += p->pSuperLib->AreaInv;
        Temp[nTemp].Score = Area + 0.01f * MAP_MAX(Delay0, Delay1);
        nTemp++;
    }
    if ( nTemp == 0 )
        return 0;
    qsort( Temp, nTemp, sizeof(Map_MctsCombo_t), Map_MctsCompareCombos );
    if ( nTemp > MAP_MCTS_MAX_COMBOS )
        nTemp = MAP_MCTS_MAX_COMBOS;
    for ( c = 0; c < nTemp; c++ )
    {
        pDecision->Combos[c] = Temp[c];
        pDecision->Combos[c].Prior = 1.0f / (float)(c + 1);
        // pDecision->Combos[c].Prior = 1.0f;
    }
    pDecision->nCombos = nTemp;
    return (pDecision->nCombos > 1);
}

static float Map_MctsDecisionPotentialScore( Map_MctsRegion_t * pRegion, Map_MctsDecision_t * pDecision )   // 综合候选丰富度、共享程度、TFI/叶子交叠和节点自身共享潜力给 decision 打分
{
    Map_Node_t * pNode0, * pNode1;
    Map_Node_t * pTfi0[MAP_MCTS_TFI_MAXNODES], * pTfi1[MAP_MCTS_TFI_MAXNODES];
    Map_Node_t * pLeaves0[2 * MAP_MCTS_MAX_PHASECANDS + 4], * pLeaves1[2 * MAP_MCTS_MAX_PHASECANDS + 4];
    float SharedRefs, PairRefs, Reconv, StructReconv, Area;
    int i, nNodes, nTfi0, nTfi1, nLeaves0, nLeaves1, PairTfiOverlap, PairLeafOverlap;
    pNode0 = pDecision->pNode;
    nNodes = Map_NodeVecReadSize( pRegion->vNodes );
    nTfi0 = 0;
    Map_MctsCollectCurrentTfi( pNode0, pTfi0, &nTfi0, MAP_MCTS_TFI_MAXNODES );
    nLeaves0 = Map_MctsCollectCurrentLeaves( pNode0, pLeaves0, 2 * MAP_MCTS_MAX_PHASECANDS + 4 );
    PairRefs = 0.0;
    PairTfiOverlap = 0;
    PairLeafOverlap = 0;
    for ( i = 0; i < nNodes; i++ )
    {
        pNode1 = Map_NodeVecReadEntry( pRegion->vNodes, i );
        if ( pNode1 == pNode0 )
            continue;
        nTfi1 = 0;
        Map_MctsCollectCurrentTfi( pNode1, pTfi1, &nTfi1, MAP_MCTS_TFI_MAXNODES );
        nLeaves1 = Map_MctsCollectCurrentLeaves( pNode1, pLeaves1, 2 * MAP_MCTS_MAX_PHASECANDS + 4 );
        PairTfiOverlap += Map_MctsCountNodeOverlap( pTfi0, nTfi0, pTfi1, nTfi1, &SharedRefs );
        PairRefs += SharedRefs;
        PairLeafOverlap += Map_MctsCountNodeOverlap( pLeaves0, nLeaves0, pLeaves1, nLeaves1, &SharedRefs );
    }
    Reconv = Map_MctsNodeReconvScore( pNode0 );
    StructReconv = Map_MctsNodeHotspotScore( pNode0 );
    Area = Map_MctsNodeArea( pNode0 );
    return 14.0f * PairTfiOverlap + 8.0f * PairLeafOverlap + PairRefs + 6.0f * pDecision->nCombos + 2.0f * pNode0->nRefs + 4.0f * Reconv + 2.0f * StructReconv + 0.2f * Area;
}

static int Map_MctsCompareDecisionScores( const void * p1, const void * p2 )
{
    Map_MctsDecision_t * pDec1 = (Map_MctsDecision_t *)p1;
    Map_MctsDecision_t * pDec2 = (Map_MctsDecision_t *)p2;
    if ( pDec1->PotScore < pDec2->PotScore )
        return 1;
    if ( pDec1->PotScore > pDec2->PotScore )
        return -1;
    return pDec2->pNode->Level - pDec1->pNode->Level;
}

static void Map_MctsBuildRegion( Map_Man_t * p, Map_Node_t * pRoot, Map_MctsRegion_t * pRegion )    // 先收集一个局部区域里的节点，再从这些节点里挑出值得做 MCTS 决策的节点，为每个节点构造 decision
{
    Map_NodeVec_t * vFanouts, * vMerges;
    Map_MctsDecision_t * pTempDecs;
    Map_Node_t * pNode;
    unsigned char * pReach;
    int i, nDecs, nTempDecs;
    memset( pRegion, 0, sizeof(Map_MctsRegion_t) );
    pRegion->pRoot = pRoot;
    pRegion->vNodes = Map_NodeVecAlloc( p->nMctsRegionMax );
    vFanouts = Map_NodeVecAlloc( 8 );
    vMerges  = Map_NodeVecAlloc( 8 );
    pReach = ABC_CALLOC( unsigned char, p->vMapObjs->nSize );
    if ( Map_MctsCollectSourceMerges( p, pRoot, pReach, vFanouts, vMerges ) > 0 )
    {
        for ( i = 0; i < Map_NodeVecReadSize(vMerges); i++ )
            Map_MctsCollectSourceRegion_rec( Map_NodeVecReadEntry(vMerges, i), pRoot, pReach, pRegion->vNodes, p->nMctsRegionMax );
        if ( Map_NodeVecReadSize(pRegion->vNodes) > 0 )
        {
            Map_MctsScoreRegion( pRegion );
            pRegion->Score += 20.0f * (float)Map_NodeVecReadSize(vMerges);
        }
    }
    ABC_FREE( pReach );
    Map_NodeVecFree( vFanouts );
    Map_NodeVecFree( vMerges );
    Map_NodeVecSortByLevel( pRegion->vNodes );
    pTempDecs = ABC_ALLOC( Map_MctsDecision_t, Abc_MaxInt(1, Map_NodeVecReadSize(pRegion->vNodes)) );
    nTempDecs = 0;
    for ( i = 0; i < Map_NodeVecReadSize(pRegion->vNodes); i++ )
    {
        pNode = Map_NodeVecReadEntry( pRegion->vNodes, i );
        if ( !Map_NodeIsAnd(pNode) || pNode->pRepr )
            continue;
        if ( Map_MctsBuildDecision( p, pNode, &pTempDecs[nTempDecs] ) )
        {
            pTempDecs[nTempDecs].PotScore = Map_MctsDecisionPotentialScore( pRegion, &pTempDecs[nTempDecs] );
            nTempDecs++;
        }
    }
    qsort( pTempDecs, nTempDecs, sizeof(Map_MctsDecision_t), Map_MctsCompareDecisionScores );
    nDecs = Abc_MinInt( nTempDecs, MAP_MCTS_MAX_DECISIONS );
    for ( i = 0; i < nDecs; i++ )
        pRegion->Decisions[i] = pTempDecs[i];
    ABC_FREE( pTempDecs );
    pRegion->nDecisions = nDecs;
}

static int Map_MctsRecomputeCover( Map_Man_t * p )    // 根据当前每个节点选中的 best cut，重新把整张映射网络的到达时间和引用计数信息计算一遍
{
    Map_Node_t * pNode;
    int i;
    Map_MappingSetPiArrivalTimes( p );
    for ( i = 0; i < p->vMapObjs->nSize; i++ )
    {
        pNode = p->vMapObjs->pArray[i];
        if ( Map_NodeIsBuf(pNode) )
        {
            pNode->tArrival[0] = Map_Regular(pNode->p1)->tArrival[ Map_IsComplement(pNode->p1) ];
            pNode->tArrival[1] = Map_Regular(pNode->p1)->tArrival[!Map_IsComplement(pNode->p1) ];
            continue;
        }
        if ( !Map_NodeIsAnd(pNode) || pNode->pRepr )
            continue;
        if ( pNode->pCutBest[0] == NULL && pNode->pCutBest[1] == NULL )
            return 0;
        if ( pNode->pCutBest[0] )
            Map_TimeCutComputeArrival( pNode, pNode->pCutBest[0], 0, MAP_FLOAT_LARGE );
        if ( pNode->pCutBest[1] )
            Map_TimeCutComputeArrival( pNode, pNode->pCutBest[1], 1, MAP_FLOAT_LARGE );
        Map_NodeTransferArrivalTimes( p, pNode );
    }
    Map_MappingSetRefs( p );
    return 1;
}

static int Map_MctsApplyChoices( Map_MctsRegion_t * pRegion, int * pChoices )    // 把外部给定的选择结果 pChoices，真正写回到 region 里对应节点上
{
    Map_Node_t * pNode;
    int i, iChoice;
    for ( i = 0; i < pRegion->nDecisions; i++ )
    {
        pNode = pRegion->Decisions[i].pNode;
        iChoice = pChoices[i] >= 0 ? pChoices[i] : 0;
        pNode->pCutBest[0] = pRegion->Decisions[i].Combos[iChoice].pCutBest[0];
        pNode->pCutBest[1] = pRegion->Decisions[i].Combos[iChoice].pCutBest[1];
        if ( pNode->pCutBest[0] == NULL && pNode->pCutBest[1] == NULL )
            return 0;
    }
    return 1;
}

static int Map_MctsCompleteChoices( Map_MctsCtx_t * pCtx, int * pChoicesIn, int * pChoicesOut )    // heuristic completion
{
    Map_MctsRegion_t * pRegion;
    Map_MctsSnapshot_t StageSnap;
    Map_Node_t * pNode;
    Map_Cut_t * pBestCut0, * pBestCut1;
    float BestArea, BestDelay, Area, Delay;
    int i, c, BestChoice;

    pRegion = pCtx->pRegion;
    memset( &StageSnap, 0, sizeof(StageSnap) );
    for ( i = 0; i < pRegion->nDecisions; i++ )
        pChoicesOut[i] = pChoicesIn[i];

    Map_MctsSnapshotRestore( pCtx->p, pCtx->pSnapshot );
    for ( i = 0; i < pRegion->nDecisions; i++ )
    {
        if ( pChoicesOut[i] < 0 )
            continue;
        pNode = pRegion->Decisions[i].pNode;
        pNode->pCutBest[0] = pRegion->Decisions[i].Combos[pChoicesOut[i]].pCutBest[0];
        pNode->pCutBest[1] = pRegion->Decisions[i].Combos[pChoicesOut[i]].pCutBest[1];
    }
    if ( !Map_MctsRecomputeCover( pCtx->p ) )
        return 0;

    for ( i = 0; i < pRegion->nDecisions; i++ )    // 开始逐个补未指定的 decision
    {
        if ( pChoicesOut[i] >= 0 )
            continue;

        Map_MctsSnapshotFree( &StageSnap );
        Map_MctsSnapshotSave( pCtx->p, &StageSnap );
        BestChoice = -1;
        BestArea   = MAP_FLOAT_LARGE;
        BestDelay  = MAP_FLOAT_LARGE;
        pBestCut0  = NULL;
        pBestCut1  = NULL;

        for ( c = 0; c < pRegion->Decisions[i].nCombos; c++ )
        {
            Map_MctsSnapshotRestore( pCtx->p, &StageSnap );
            pNode = pRegion->Decisions[i].pNode;
            pNode->pCutBest[0] = pRegion->Decisions[i].Combos[c].pCutBest[0];
            pNode->pCutBest[1] = pRegion->Decisions[i].Combos[c].pCutBest[1];
            if ( !Map_MctsRecomputeCover( pCtx->p ) )
                continue;

            Delay = Map_TimeComputeArrivalMax( pCtx->p );
            Area  = Map_MappingGetArea( pCtx->p );
            if ( Delay > pCtx->DelayBound + pCtx->p->fEpsilon )
                continue;

            if ( BestChoice == -1 || Area < BestArea - pCtx->p->fEpsilon ||
                (Area < BestArea + pCtx->p->fEpsilon && Delay < BestDelay - pCtx->p->fEpsilon) )    // 第一优先级：更小 area；第二优先级：更小 delay
            {
                BestChoice = c;
                BestArea   = Area;
                BestDelay  = Delay;
                pBestCut0  = pNode->pCutBest[0];
                pBestCut1  = pNode->pCutBest[1];
            }
            else if ( Area < BestArea + pCtx->p->fEpsilon && Delay < BestDelay + pCtx->p->fEpsilon && BestChoice >= 0 )    // 第二层 tie-break
            {
                int fReplace = 0;
                if ( pBestCut0 && pNode->pCutBest[0] )
                    fReplace |= Map_MatchCompare( pCtx->p, &pBestCut0->M[0], &pNode->pCutBest[0]->M[0], 1 );
                if ( pBestCut1 && pNode->pCutBest[1] )
                    fReplace |= Map_MatchCompare( pCtx->p, &pBestCut1->M[1], &pNode->pCutBest[1]->M[1], 1 );
                if ( fReplace )
                {
                    BestChoice = c;
                    BestArea   = Area;
                    BestDelay  = Delay;
                    pBestCut0  = pNode->pCutBest[0];
                    pBestCut1  = pNode->pCutBest[1];
                }
            }
        }

        if ( BestChoice < 0 )   // 如果所有 combo 都不合法，则失败
        {
            Map_MctsSnapshotFree( &StageSnap );
            return 0;
        }

        pChoicesOut[i] = BestChoice;    // 找到最优 combo 后，把它正式写入输出和当前状态
        Map_MctsSnapshotRestore( pCtx->p, &StageSnap );
        pNode = pRegion->Decisions[i].pNode;
        pNode->pCutBest[0] = pRegion->Decisions[i].Combos[BestChoice].pCutBest[0];
        pNode->pCutBest[1] = pRegion->Decisions[i].Combos[BestChoice].pCutBest[1];
        if ( !Map_MctsRecomputeCover( pCtx->p ) )
        {
            Map_MctsSnapshotFree( &StageSnap );
            return 0;
        }
    }
    Map_MctsSnapshotFree( &StageSnap );
    return 1;
}

static float Map_MctsEvaluateChoices( Map_MctsCtx_t * pCtx, int * pChoices )
{
    Map_MctsRegion_t * pRegion;
    float Area, Delay, Value;
    int Choices[MAP_MCTS_MAX_DECISIONS];
    int fValid;
    int i;

    pRegion = pCtx->pRegion;
    // Region search runs on a private manager per worker, so evaluation does not
    // need a global critical section. Keeping it serialized defeats parallel MCTS.
    fValid = Map_MctsCompleteChoices( pCtx, pChoices, Choices );    // 先把“候选动作”变成“一个完整可执行的方案”
    if ( fValid )
    {
        Delay = Map_TimeComputeArrivalMax( pCtx->p );
        Area  = Map_MappingGetArea( pCtx->p );
    }
    Map_MctsSnapshotRestore( pCtx->p, pCtx->pSnapshot );
    if ( !fValid )
        return -10000.0;
    if ( Delay > pCtx->DelayBound + pCtx->p->fEpsilon )
        Value = -1000.0f - 10.0f * (Delay - pCtx->DelayBound);
    else
        Value = pCtx->AreaBase - Area;
    if ( Value > pCtx->BestValue )
    {
        pCtx->BestValue = Value;
        pCtx->BestArea  = Area;
        pCtx->BestDelay = Delay;
        for ( i = 0; i < pRegion->nDecisions; i++ )
            pCtx->BestChoices[i] = Choices[i];
        if ( pCtx->p->fVerbose )
        {
            Map_MctsTrace( "  region %d sim %d best update: value=%8.2f area=%11.1f delay=%8.2f",
                pCtx->iRegion + 1, pCtx->iSim + 1, Value, Area, Delay );
            Map_MctsPrintChoices( pRegion, Choices );
        }
    }
    return Value;
}

static float Map_MctsSearch( Map_MctsCtx_t * pCtx, Map_MctsTree_t * pTree, int * pChoices )    // 选择 → 扩展/评估 → 回传
{
    Map_MctsDecision_t * pDecision;
    Map_MctsTree_t * pChild;
    float BestScore, Score, Q;
    float Value;
    int BestChild, i;

    if ( pTree->Depth == pCtx->pRegion->nDecisions )           // 终止1：深度已经等于 decision 数
    {
        Value = Map_MctsEvaluateChoices( pCtx, pChoices );
        pTree->nVisits++;
        pTree->ValueSum += Value;
        return Value;
    }
    if ( pTree->nVisits == 0 )                                // 终止2：当前节点从未访问过
    {
        Value = Map_MctsEvaluateChoices( pCtx, pChoices );
        pTree->nVisits++;
        pTree->ValueSum += Value;
        return Value;
    }

    if ( pTree->nChildren == 0 )                              // 未访问：初始化 children 和先验
    {
        pDecision = &pCtx->pRegion->Decisions[pTree->Depth];
        pTree->nChildren = pDecision->nCombos;
        for ( i = 0; i < pTree->nChildren; i++ )
            pTree->Priors[i] = pDecision->Combos[i].Prior;
    }

    BestChild = 0;
    BestScore = -MAP_FLOAT_LARGE;
    for ( i = 0; i < pTree->nChildren; i++ )                  // selection
    {
        pChild = pTree->pChildren[i];
        if ( pChild == NULL )
            Q = 0.0;
        else
            Q = pChild->ValueSum / (float)Abc_MaxInt( pChild->nVisits, 1 );
        Score = Q + 1.5f * pTree->Priors[i] * sqrt((double)Abc_MaxInt(pTree->nVisits, 1)) / (1.0f + (float)(pChild ? pChild->nVisits : 0));
        if ( Score > BestScore )
            BestScore = Score, BestChild = i;
    }

    if ( pTree->pChildren[BestChild] == NULL )
    {
        pTree->pChildren[BestChild] = ABC_CALLOC( Map_MctsTree_t, 1 );
        pTree->pChildren[BestChild]->Depth = pTree->Depth + 1;
    }
    pChoices[pTree->Depth] = BestChild;
    Value = Map_MctsSearch( pCtx, pTree->pChildren[BestChild], pChoices );
    pChoices[pTree->Depth] = -1;
    pTree->nVisits++;                                        // backup
    pTree->ValueSum += Value;
    return Value;
}

static int Map_MctsSearchOneRegion( Map_Man_t * p, Map_Node_t * pRoot, float AreaBase, float DelayBase, int iRegion, int nRegions, Map_MctsResult_t * pRes )
{
    Map_MctsSnapshot_t Snapshot;
    Map_MctsRegion_t Region;
    Map_MctsCtx_t Ctx;
    Map_MctsTree_t * pTree;
    int Choices[MAP_MCTS_MAX_DECISIONS];
    int k;

    memset( pRes, 0, sizeof(Map_MctsResult_t) );
    memset( &Snapshot, 0, sizeof(Snapshot) );
    memset( &Region, 0, sizeof(Region) );
    Map_MctsSnapshotSave( p, &Snapshot );
    Map_MctsBuildRegion( p, pRoot, &Region );   // 构建当前 region
    if ( Region.nDecisions == 0 )
    {
        Map_NodeVecFree( Region.vNodes );
        Map_MctsSnapshotFree( &Snapshot );
        return 0;
    }

    memset( &Ctx, 0, sizeof(Ctx) );
    Ctx.p         = p;
    Ctx.pRegion   = &Region;
    Ctx.pSnapshot = &Snapshot;
    Ctx.pRoot     = pRoot;
    Ctx.AreaBase  = AreaBase;
    Ctx.DelayBound = DelayBase;
    Ctx.BestValue = -MAP_FLOAT_LARGE;
    Ctx.BestArea  = AreaBase;
    Ctx.BestDelay = DelayBase;
    Ctx.iRegion   = iRegion;
    Ctx.nRegions  = nRegions;
    for ( k = 0; k < MAP_MCTS_MAX_DECISIONS; k++ )
        Ctx.BestChoices[k] = 0, Choices[k] = -1;

    if ( p->fVerbose )
        Map_MctsPrintRegion( &Ctx );
    pTree = ABC_CALLOC( Map_MctsTree_t, 1 );        // 创建 MCTS 树根
    for ( k = 0; k < p->nMctsSimNum; k++ )          // 做多轮 simulation
    {
        Ctx.iSim = k;
        Map_MctsSearch( &Ctx, pTree, Choices );
    }
    Map_MctsTreeFree( pTree );

    if ( Ctx.BestValue > p->fEpsilon )              // 如果找到正收益改进，就写入结果
    {
        pRes->fFound = 1;
        pRes->BestValue = Ctx.BestValue;
        pRes->BestArea = Ctx.BestArea;
        pRes->BestDelay = Ctx.BestDelay;
        for ( k = 0; k < MAP_MCTS_MAX_DECISIONS; k++ )
            pRes->BestChoices[k] = Ctx.BestChoices[k];
        if ( p->fVerbose )
        {
            Map_MctsTrace( "  region %d candidate: value=%8.2f area=%11.1f delay=%8.2f",
                iRegion + 1, Ctx.BestValue, Ctx.BestArea, Ctx.BestDelay );
            Map_MctsPrintChoices( &Region, Ctx.BestChoices );
        }
    }
    else if ( p->fVerbose )
        Map_MctsTrace( "  region %d skipped: no improving assignment found\n", iRegion + 1 );

    Map_NodeVecFree( Region.vNodes );
    Map_MctsSnapshotFree( &Snapshot );
    return 1;
}

static void Map_MctsTreeFree( Map_MctsTree_t * pTree )
{
    int i;
    if ( pTree == NULL )
        return;
    for ( i = 0; i < pTree->nChildren; i++ )
        Map_MctsTreeFree( pTree->pChildren[i] );
    ABC_FREE( pTree );
}

int Map_MappingMctsRefine( Map_Man_t * p )     // 枚举所有 hotspot region，并行/串行执行
{
    Map_MctsSnapshot_t Snapshot;
    Map_MctsRegion_t * pRegions;
    Map_MctsResult_t * pResults;
    float AreaBase, DelayBase, AreaNew;
    float DelayNew;
    int i, nHot, nAccepted, nRegions, nRegionCands, nRegionCandsMax;
    abctime clk;

    if ( !p->fMctsEnable || p->nMctsRegionNum <= 0 || p->nMctsSimNum <= 0 )
        return 0;

    clk = Abc_Clock();
    if ( p->fVerbose )
        Map_MctsTraceOpen();
    Map_MctsSnapshotSave( p, &Snapshot );
    AreaBase  = Map_MappingGetArea( p );
    DelayBase = (p->DelayTarget >= 0.0 && p->DelayTarget < ABC_INFINITY) ? p->DelayTarget : Map_TimeComputeArrivalMax( p );
    nHot = 0;
    nRegionCands = 0;
    nRegionCandsMax = Abc_MinInt( p->vMapObjs->nSize, Abc_MaxInt( p->nMctsRegionNum, 16 * p->nMctsRegionNum ) );
    pRegions = ABC_CALLOC( Map_MctsRegion_t, nRegionCandsMax );
    pResults = NULL;

    if ( p->fVerbose )
    {
        Mio_Library_t * pLib = p->pSuperLib ? p->pSuperLib->pGenlib : NULL;
        Map_MctsTrace( "MCTS context: file=%s lib=%s lib_file=%s super=%s\n",
            Map_MctsLogStr(p->pFileName),
            Map_MctsLogStr(pLib ? Mio_LibraryReadName(pLib) : NULL),
            Map_MctsLogStr(pLib ? Mio_LibraryReadFileName(pLib) : NULL),
            Map_MctsLogStr(p->pSuperLib ? p->pSuperLib->pName : NULL) );
        Map_MctsTrace( "MCTS start: scanning reconvergent sources region_limit=%d candidate_limit=%d sims_per_region=%d max_region_nodes=%d base_area=%11.1f delay_bound=%8.2f log_file=%s\n",
            p->nMctsRegionNum, nRegionCandsMax, p->nMctsSimNum, p->nMctsRegionMax, AreaBase, DelayBase, MAP_MCTS_LOG_FILE );
        Map_MctsTrace( "MCTS region search mode: %s threads=%d\n",
            p->fMctsParallel ?
#ifdef ABC_USE_OPENMP
            "parallel"
#else
            "serial"
#endif
            : "serial",
            p->fMctsParallel ?
#ifdef ABC_USE_OPENMP
            omp_get_max_threads() : 1 );
#else
            1 : 1 );
#endif
    }

    nAccepted = 0;
    for ( i = 0; i < p->vMapObjs->nSize; i++ )
    {
        Map_MctsRegion_t Region;
        int k, iWorst;
        Map_Node_t * pNode = p->vMapObjs->pArray[i];
        if ( !Map_NodeIsAnd(pNode) || pNode->pRepr )
            continue;
        if ( !Map_MctsIsReconvergentSource( p, pNode ) )
            continue;
        nHot++;
        memset( &Region, 0, sizeof(Region) );
        Map_MctsBuildRegion( p, pNode, &Region );
        if ( Map_NodeVecReadSize(Region.vNodes) == 0 || Region.Score <= MAP_FLOAT_SMALL )
        {
            Map_NodeVecFree( Region.vNodes );
            continue;
        }
        if ( nRegionCands < nRegionCandsMax )
        {
            pRegions[nRegionCands++] = Region;
            continue;
        }
        iWorst = 0;
        for ( k = 1; k < nRegionCands; k++ )
            if ( Map_MctsCompareRegions( &pRegions[k], &pRegions[iWorst] ) > 0 )
                iWorst = k;
        if ( Map_MctsCompareRegions( &Region, &pRegions[iWorst] ) < 0 )
        {
            Map_NodeVecFree( pRegions[iWorst].vNodes );
            pRegions[iWorst] = Region;
        }
        else
            Map_NodeVecFree( Region.vNodes );
    }
    nRegions = Abc_MinInt( nRegionCands, p->nMctsRegionNum );
    pResults = ABC_CALLOC( Map_MctsResult_t, nRegions );
    if ( p->fVerbose )
        Map_MctsTrace( "MCTS sources: reconvergent=%d retained_regions=%d\n", nHot, nRegionCands );
    qsort( pRegions, nRegionCands, sizeof(Map_MctsRegion_t), Map_MctsCompareRegions );

    if ( !p->fMctsParallel )                                    // 对每个 region 跑搜索：串行或并行
    {
        for ( i = 0; i < nRegions; i++ )
        {
            if ( pRegions[i].nDecisions == 0 )
                continue;
            Map_MctsSearchOneRegion( p, pRegions[i].pRoot, AreaBase, DelayBase, i, nRegions, &pResults[i] );
        }
    }
    else
    {
#ifdef ABC_USE_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
        for ( i = 0; i < nRegions; i++ )
        {
            Map_Man_t * pLocal;
            if ( pRegions[i].nDecisions == 0 )
                continue;
            pLocal = Map_MctsManDup( p );
            Map_MctsSearchOneRegion( pLocal, pLocal->vMapObjs->pArray[pRegions[i].pRoot->Num], AreaBase, DelayBase, i, nRegions, &pResults[i] );
            Map_MctsManFree( pLocal );
        }
    }

    for ( i = 0; i < nRegions; i++ )        // 搜索结束后逐个 replay 验证
    {
        if ( pResults[i].fFound )
        {
            Map_MctsSnapshotRestore( p, &Snapshot );
            if ( Map_MctsApplyChoices( &pRegions[i], pResults[i].BestChoices ) && Map_MctsRecomputeCover(p) )
            {
                AreaNew = Map_MappingGetArea( p );
                DelayNew = Map_TimeComputeArrivalMax( p );
                if ( DelayNew <= DelayBase + p->fEpsilon && AreaNew + p->fEpsilon < AreaBase )
                {
                    Map_TimeComputeRequiredGlobal( p );
                    Snapshot.AreaFinal = p->AreaFinal = AreaNew;
                    Snapshot.RequiredGlo = p->fRequiredGlo;
                    Map_MctsSnapshotFree( &Snapshot );
                    Map_MctsSnapshotSave( p, &Snapshot );
                    if ( p->fVerbose )
                        Map_MctsTrace( "  region %d accepted: area %11.1f -> %11.1f  delay=%8.2f\n",
                            i + 1, AreaBase, AreaNew, DelayNew );
                    AreaBase = AreaNew;
                    nAccepted++;
                }
                else
                {
                    if ( p->fVerbose )
                        Map_MctsTrace( "  region %d rejected after replay: area=%11.1f delay=%8.2f (base_area=%11.1f bound=%8.2f)\n",
                            i + 1, AreaNew, DelayNew, AreaBase, DelayBase );
                    Map_MctsSnapshotRestore( p, &Snapshot );
                }
            }
            else if ( p->fVerbose )
            {
                if ( p->fVerbose )
                    Map_MctsTrace( "  region %d rejected: replay failed\n", i + 1 );
            }
        }
        Map_MctsSnapshotRestore( p, &Snapshot );
        Map_NodeVecFree( pRegions[i].vNodes );
    }
    for ( ; i < nRegionCands; i++ )
        Map_NodeVecFree( pRegions[i].vNodes );
    ABC_FREE( pRegions );
    ABC_FREE( pResults );
    Map_MctsSnapshotRestore( p, &Snapshot );
    Map_MctsSnapshotFree( &Snapshot );

    if ( nAccepted > 0 )
    {
        Map_MappingSetRefs( p );
        p->AreaFinal = Map_MappingGetArea( p );
        Map_TimeComputeRequiredGlobal( p );
    }
    if ( p->fVerbose )
        Map_MctsTrace( "MCTS summary: accepted=%d runtime=%9.2f sec final_area=%11.1f final_delay=%8.2f\n",
            nAccepted, 1.0*((double)(Abc_Clock() - clk))/((double)CLOCKS_PER_SEC), p->AreaFinal, p->fRequiredGlo );
    Map_MctsTraceClose();
    return nAccepted;
}

ABC_NAMESPACE_IMPL_END
