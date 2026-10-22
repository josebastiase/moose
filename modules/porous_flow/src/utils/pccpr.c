/*
  PCCPR: two-stage Constrained Pressure Residual preconditioner with
  quasi-IMPES decoupling.

      stage 1:  y1 = P * Ksp_p( R r )          Ap = R A P   (AMG)
      stage 2:  y  = y1 + Ksp_g( r - A y1 )                 (ILU / bjacobi)

  R has one row per node.  The weights in that row decide how the node's
  equations are combined into a single "pressure equation":
      qimpes : w_i = D_ii^{-T} e_p   -> zero coupling to the node's own
                                        non-pressure unknowns (default)
      sum    : w_i = (1,...,1)       -> sum of component mass balances
      raw    : w_i = e_p             -> just the pressure-variable row
  P injects a pressure correction into the pressure dof of each node.

  Author: Jose (SubTerra / PorousFlow work), with Claude.  BSD-2-Clause.
*/
#include <petsc/private/pcimpl.h>
#include <petscksp.h>
#include "pccpr.h"

#if PETSC_VERSION_LT(3, 23, 0)
  #define CPR_OPTITEMS PetscOptionItems *
#else
  #define CPR_OPTITEMS PetscOptionItems
#endif

#if PETSC_VERSION_GE(3, 20, 0)
  #define CPR_NEST(ksp, pc) PetscCall(KSPSetNestLevel(ksp, (pc)->kspnestlevel))
#else
  #define CPR_NEST(ksp, pc) (void)0
#endif

typedef enum {
  PC_CPR_WEIGHT_QIMPES,
  PC_CPR_WEIGHT_SUM,
  PC_CPR_WEIGHT_RAW
} PCCPRWeightType;
static const char *const PCCPRWeightTypes[] = {"QIMPES", "SUM", "RAW", "PCCPRWeightType", "PC_CPR_WEIGHT_", NULL};

typedef enum {
  PC_CPR_MULTIPLICATIVE,
  PC_CPR_ADDITIVE
} PCCPRCompositeType;
static const char *const PCCPRCompositeTypes[] = {"MULTIPLICATIVE", "ADDITIVE", "PCCPRCompositeType", "PC_CPR_", NULL};

typedef struct {
  PetscInt           bs;     /* dofs per node; 0 = take from matrix */
  PetscInt           pfield; /* pressure position within a node     */
  PCCPRWeightType    wtype;
  PCCPRCompositeType ctype;

  /* layout: nloc nodes, dofs[k*bs + f] = global index of field f at node k */
  PetscInt  nloc, lbs, nstart;
  PetscInt *dofs;
  IS        layout_is; /* user-supplied layout this setup was built from */

  Mat       R, P, Ap;
  KSP       kspP, kspG;
  Vec       rp, xp, r2, z;
  PetscInt  nfallback; /* nodes where the decoupling block was singular */
} PC_CPR;

/* ---------------------------------------------------------------- helpers */

/* solve M x = b in place (M n x n row-major), partial pivoting */
static PetscBool CPRSmallSolve(PetscInt n, PetscScalar *M, PetscScalar *b)
{
  for (PetscInt k = 0; k < n; k++) {
    PetscInt  piv = k;
    PetscReal big = PetscAbsScalar(M[k * n + k]);
    for (PetscInt i = k + 1; i < n; i++)
      if (PetscAbsScalar(M[i * n + k]) > big) {
        big = PetscAbsScalar(M[i * n + k]);
        piv = i;
      }
    if (big == 0.0) return PETSC_FALSE;
    if (piv != k) {
      for (PetscInt j = 0; j < n; j++) {
        PetscScalar t = M[k * n + j];
        M[k * n + j]  = M[piv * n + j];
        M[piv * n + j] = t;
      }
      PetscScalar t = b[k];
      b[k]          = b[piv];
      b[piv]        = t;
    }
    for (PetscInt i = k + 1; i < n; i++) {
      PetscScalar f = M[i * n + k] / M[k * n + k];
      for (PetscInt j = k; j < n; j++) M[i * n + j] -= f * M[k * n + j];
      b[i] -= f * b[k];
    }
  }
  for (PetscInt k = n - 1; k >= 0; k--) {
    for (PetscInt j = k + 1; j < n; j++) b[k] -= M[k * n + j] * b[j];
    b[k] /= M[k * n + k];
  }
  for (PetscInt k = 0; k < n; k++)
    if (PetscIsInfOrNanScalar(b[k])) return PETSC_FALSE;
  return PETSC_TRUE;
}

static PetscErrorCode CPRCreateKSPs(PC pc)
{
  PC_CPR     *cpr = (PC_CPR *)pc->data;
  const char *prefix;
  PC          sub;

  PetscFunctionBegin;
  PetscCall(PCGetOptionsPrefix(pc, &prefix));
  if (!cpr->kspP) {
    PetscCall(KSPCreate(PetscObjectComm((PetscObject)pc), &cpr->kspP));
    CPR_NEST(cpr->kspP, pc);
    PetscCall(PetscObjectIncrementTabLevel((PetscObject)cpr->kspP, (PetscObject)pc, 1));
    PetscCall(KSPSetOptionsPrefix(cpr->kspP, prefix));
    PetscCall(KSPAppendOptionsPrefix(cpr->kspP, "pc_cpr_pressure_"));
    PetscCall(KSPSetType(cpr->kspP, KSPPREONLY));
    PetscCall(KSPGetPC(cpr->kspP, &sub));
#if defined(PETSC_HAVE_HYPRE)
    PetscCall(PCSetType(sub, PCHYPRE));
    PetscCall(PCHYPRESetType(sub, "boomeramg"));
#else
    PetscCall(PCSetType(sub, PCGAMG));
#endif
  }
  if (!cpr->kspG) {
    PetscCall(KSPCreate(PetscObjectComm((PetscObject)pc), &cpr->kspG));
    CPR_NEST(cpr->kspG, pc);
    PetscCall(PetscObjectIncrementTabLevel((PetscObject)cpr->kspG, (PetscObject)pc, 1));
    PetscCall(KSPSetOptionsPrefix(cpr->kspG, prefix));
    PetscCall(KSPAppendOptionsPrefix(cpr->kspG, "pc_cpr_global_"));
    PetscCall(KSPSetType(cpr->kspG, KSPPREONLY));
    PetscCall(KSPGetPC(cpr->kspG, &sub));
    PetscCall(PCSetType(sub, PCBJACOBI));
    {
      /* Default to a NONZERO shift in the stage-2 sub-block factorisations.
         PorousFlow transport/tracer rows can have (near-)zero pivots; without
         the shift ILU(0) fails with SUBPC_ERROR.  User options still win. */
      char      opt[PETSC_MAX_PATH_LEN];
      PetscBool set;
      PetscCall(PetscSNPrintf(opt, sizeof(opt), "-%spc_cpr_global_sub_pc_factor_shift_type", prefix ? prefix : ""));
      PetscCall(PetscOptionsHasName(((PetscObject)pc)->options, NULL, opt, &set));
      if (!set) PetscCall(PetscOptionsSetValue(((PetscObject)pc)->options, opt, "NONZERO"));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* decide which global dofs belong to which node */
static PetscErrorCode CPRBuildLayout(PC pc, IS user)
{
  PC_CPR  *cpr = (PC_CPR *)pc->data;
  PetscInt rstart, rend, mloc, bs, nend;

  PetscFunctionBegin;
  PetscCall(PetscFree(cpr->dofs));
  PetscCall(MatGetOwnershipRange(pc->pmat, &rstart, &rend));
  mloc = rend - rstart;
  if (user) {
    const PetscInt *idx;
    PetscInt        n;
    PetscCall(ISGetBlockSize(user, &bs));
    PetscCall(ISGetLocalSize(user, &n));
    PetscCall(ISGetIndices(user, &idx));
    cpr->nloc = n / bs;
    PetscCall(PetscMalloc1(n, &cpr->dofs));
    for (PetscInt i = 0; i < n; i++) {
      PetscCheck(idx[i] >= rstart && idx[i] < rend, PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE,
                 "PCCPR: node dof %" PetscInt_FMT " is not owned by this rank [%" PetscInt_FMT ",%" PetscInt_FMT ")", idx[i], rstart, rend);
      cpr->dofs[i] = idx[i];
    }
    PetscCall(ISRestoreIndices(user, &idx));
  } else {
    bs = cpr->bs;
    if (bs <= 0) PetscCall(MatGetBlockSize(pc->pmat, &bs));
    PetscCheck(bs >= 1, PetscObjectComm((PetscObject)pc), PETSC_ERR_ARG_WRONG, "PCCPR: set -pc_cpr_block_size");
    PetscCheck(mloc % bs == 0, PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ,
               "PCCPR: local rows %" PetscInt_FMT " not divisible by block size %" PetscInt_FMT " (dofs not interleaved per node?)", mloc, bs);
    cpr->nloc = mloc / bs;
    PetscCall(PetscMalloc1(mloc, &cpr->dofs));
    for (PetscInt i = 0; i < mloc; i++) cpr->dofs[i] = rstart + i;
  }
  PetscCheck(cpr->pfield >= 0 && cpr->pfield < bs, PetscObjectComm((PetscObject)pc), PETSC_ERR_ARG_OUTOFRANGE,
             "PCCPR: pressure field %" PetscInt_FMT " outside block size %" PetscInt_FMT, cpr->pfield, bs);
  cpr->lbs = bs;
  PetscCallMPI(MPI_Scan(&cpr->nloc, &nend, 1, MPIU_INT, MPI_SUM, PetscObjectComm((PetscObject)pc)));
  cpr->nstart = nend - cpr->nloc;
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* fill R with the current weights */
static PetscErrorCode CPRComputeWeights(PC pc)
{
  PC_CPR      *cpr = (PC_CPR *)pc->data;
  PetscInt     bs  = cpr->lbs, nfb = 0;
  PetscScalar *D, *Mt, *w;

  PetscFunctionBegin;
  PetscCall(PetscMalloc3(bs * bs, &D, bs * bs, &Mt, bs, &w));
  for (PetscInt k = 0; k < cpr->nloc; k++) {
    const PetscInt *nd  = cpr->dofs + k * bs;
    PetscInt        row = cpr->nstart + k;

    for (PetscInt f = 0; f < bs; f++) w[f] = 0.0;
    switch (cpr->wtype) {
    case PC_CPR_WEIGHT_RAW:
      w[cpr->pfield] = 1.0;
      break;
    case PC_CPR_WEIGHT_SUM:
      for (PetscInt f = 0; f < bs; f++) w[f] = 1.0;
      break;
    case PC_CPR_WEIGHT_QIMPES: {
      PetscReal wmax = 0.0;
      PetscCall(MatGetValues(pc->pmat, bs, nd, bs, nd, D)); /* D[a*bs+b] = A(row a, col b) */
      for (PetscInt a = 0; a < bs; a++)
        for (PetscInt b = 0; b < bs; b++) Mt[a * bs + b] = D[b * bs + a]; /* D^T */
      w[cpr->pfield] = 1.0;
      if (CPRSmallSolve(bs, Mt, w)) {
        for (PetscInt f = 0; f < bs; f++) wmax = PetscMax(wmax, PetscAbsScalar(w[f]));
      }
      if (wmax > 0.0) {
        for (PetscInt f = 0; f < bs; f++) w[f] /= wmax;
      } else { /* singular local block: fall back to the raw pressure row */
        for (PetscInt f = 0; f < bs; f++) w[f] = 0.0;
        w[cpr->pfield] = 1.0;
        nfb++;
      }
    } break;
    }
    PetscCall(MatSetValues(cpr->R, 1, &row, bs, nd, w, INSERT_VALUES));
  }
  PetscCall(MatAssemblyBegin(cpr->R, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(cpr->R, MAT_FINAL_ASSEMBLY));
  PetscCall(PetscFree3(D, Mt, w));
  PetscCallMPI(MPIU_Allreduce(&nfb, &cpr->nfallback, 1, MPIU_INT, MPI_SUM, PetscObjectComm((PetscObject)pc)));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* ------------------------------------------------------------- PC methods */

static PetscErrorCode PCReset_CPR(PC pc)
{
  PC_CPR *cpr = (PC_CPR *)pc->data;

  PetscFunctionBegin;
  PetscCall(MatDestroy(&cpr->R));
  PetscCall(MatDestroy(&cpr->P));
  PetscCall(MatDestroy(&cpr->Ap));
  PetscCall(VecDestroy(&cpr->rp));
  PetscCall(VecDestroy(&cpr->xp));
  PetscCall(VecDestroy(&cpr->r2));
  PetscCall(VecDestroy(&cpr->z));
  PetscCall(ISDestroy(&cpr->layout_is));
  PetscCall(PetscFree(cpr->dofs));
  if (cpr->kspP) PetscCall(KSPReset(cpr->kspP));
  if (cpr->kspG) PetscCall(KSPReset(cpr->kspG));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PCSetUp_CPR(PC pc)
{
  PC_CPR   *cpr = (PC_CPR *)pc->data;
  IS        user = NULL;
  PetscBool rebuild, isaij;

  PetscFunctionBegin;
  PetscCall(PetscObjectBaseTypeCompareAny((PetscObject)pc->pmat, &isaij, MATSEQAIJ, MATMPIAIJ, MATAIJ, ""));
  PetscCheck(isaij, PetscObjectComm((PetscObject)pc), PETSC_ERR_SUP, "PCCPR needs an assembled AIJ preconditioning matrix");
  PetscCall(PetscObjectQuery((PetscObject)pc, "PCCPR_NodeDofs", (PetscObject *)&user));

  rebuild = (PetscBool)(!cpr->R || pc->flag == DIFFERENT_NONZERO_PATTERN || user != cpr->layout_is);
  if (rebuild) {
    PetscInt mloc, nloc, rstart;

    PetscCall(MatDestroy(&cpr->R));
    PetscCall(MatDestroy(&cpr->P));
    PetscCall(MatDestroy(&cpr->Ap));
    PetscCall(VecDestroy(&cpr->rp));
    PetscCall(VecDestroy(&cpr->xp));
    PetscCall(VecDestroy(&cpr->r2));
    PetscCall(VecDestroy(&cpr->z));
    if (user) PetscCall(PetscObjectReference((PetscObject)user));
    PetscCall(ISDestroy(&cpr->layout_is));
    cpr->layout_is = user;

    PetscCall(CPRBuildLayout(pc, user));
    PetscCall(MatGetLocalSize(pc->pmat, &mloc, NULL));
    PetscCall(MatGetOwnershipRange(pc->pmat, &rstart, NULL));
    nloc = cpr->nloc;

    /* R: nodes x dofs */
    PetscCall(MatCreate(PetscObjectComm((PetscObject)pc), &cpr->R));
    PetscCall(MatSetSizes(cpr->R, nloc, mloc, PETSC_DETERMINE, PETSC_DETERMINE));
    PetscCall(MatSetType(cpr->R, MATAIJ));
    PetscCall(MatSeqAIJSetPreallocation(cpr->R, cpr->lbs, NULL));
    PetscCall(MatMPIAIJSetPreallocation(cpr->R, cpr->lbs, NULL, 0, NULL));

    /* P: dofs x nodes, a single 1 on the pressure dof of each node */
    PetscCall(MatCreate(PetscObjectComm((PetscObject)pc), &cpr->P));
    PetscCall(MatSetSizes(cpr->P, mloc, nloc, PETSC_DETERMINE, PETSC_DETERMINE));
    PetscCall(MatSetType(cpr->P, MATAIJ));
    PetscCall(MatSeqAIJSetPreallocation(cpr->P, 1, NULL));
    PetscCall(MatMPIAIJSetPreallocation(cpr->P, 1, NULL, 0, NULL));
    for (PetscInt k = 0; k < nloc; k++) {
      PetscInt    i   = cpr->dofs[k * cpr->lbs + cpr->pfield], j = cpr->nstart + k;
      PetscScalar one = 1.0;
      PetscCall(MatSetValues(cpr->P, 1, &i, 1, &j, &one, INSERT_VALUES));
    }
    PetscCall(MatAssemblyBegin(cpr->P, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(cpr->P, MAT_FINAL_ASSEMBLY));

    PetscCall(CPRComputeWeights(pc));
    PetscCall(MatMatMatMult(cpr->R, pc->pmat, cpr->P, MAT_INITIAL_MATRIX, PETSC_DETERMINE, &cpr->Ap));

    PetscCall(MatCreateVecs(cpr->Ap, &cpr->xp, &cpr->rp));
    PetscCall(MatCreateVecs(pc->pmat, &cpr->z, &cpr->r2));
  } else {
    /* same pattern: new weights, numeric triple product only */
    PetscCall(CPRComputeWeights(pc));
    PetscCall(MatMatMatMult(cpr->R, pc->pmat, cpr->P, MAT_REUSE_MATRIX, PETSC_DETERMINE, &cpr->Ap));
  }
  if (cpr->nfallback) PetscCall(PetscInfo(pc, "PCCPR: %" PetscInt_FMT " nodes used raw weights (singular diagonal block)\n", cpr->nfallback));

  PetscCall(CPRCreateKSPs(pc));
  PetscCall(KSPSetOperators(cpr->kspP, cpr->Ap, cpr->Ap));
  PetscCall(KSPSetOperators(cpr->kspG, pc->mat, pc->pmat));
  if (!pc->setupcalled || rebuild) {
    PetscCall(KSPSetFromOptions(cpr->kspP));
    PetscCall(KSPSetFromOptions(cpr->kspG));
  }
  PetscCall(KSPSetUp(cpr->kspP));
  PetscCall(KSPSetUp(cpr->kspG));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PCApply_CPR(PC pc, Vec x, Vec y)
{
  PC_CPR *cpr = (PC_CPR *)pc->data;

  PetscFunctionBegin;
  /* stage 1: decoupled pressure solve */
  PetscCall(MatMult(cpr->R, x, cpr->rp));
  PetscCall(KSPSolve(cpr->kspP, cpr->rp, cpr->xp));
  PetscCall(KSPCheckSolve(cpr->kspP, pc, cpr->xp));
  PetscCall(MatMult(cpr->P, cpr->xp, y));
  /* stage 2: global smoother on the remaining residual */
  if (cpr->ctype == PC_CPR_MULTIPLICATIVE) {
    PetscCall(MatMult(pc->pmat, y, cpr->r2));
    PetscCall(VecAYPX(cpr->r2, -1.0, x));
  } else {
    PetscCall(VecCopy(x, cpr->r2));
  }
  PetscCall(KSPSolve(cpr->kspG, cpr->r2, cpr->z));
  PetscCall(KSPCheckSolve(cpr->kspG, pc, cpr->z));
  PetscCall(VecAXPY(y, 1.0, cpr->z));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PCSetFromOptions_CPR(PC pc, CPR_OPTITEMS PetscOptionsObject)
{
  PC_CPR *cpr = (PC_CPR *)pc->data;

  PetscFunctionBegin;
  PetscOptionsHeadBegin(PetscOptionsObject, "CPR options");
  PetscCall(PetscOptionsInt("-pc_cpr_block_size", "Dofs per node (0: matrix block size)", "PCCPRSetBlockSize", cpr->bs, &cpr->bs, NULL));
  PetscCall(PetscOptionsInt("-pc_cpr_pressure_field", "Position of pressure within a node", "PCCPRSetPressureField", cpr->pfield, &cpr->pfield, NULL));
  PetscCall(PetscOptionsEnum("-pc_cpr_weights", "Decoupling weights", "", PCCPRWeightTypes, (PetscEnum)cpr->wtype, (PetscEnum *)&cpr->wtype, NULL));
  PetscCall(PetscOptionsEnum("-pc_cpr_type", "How the two stages are combined", "", PCCPRCompositeTypes, (PetscEnum)cpr->ctype, (PetscEnum *)&cpr->ctype, NULL));
  PetscOptionsHeadEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PCView_CPR(PC pc, PetscViewer viewer)
{
  PC_CPR   *cpr = (PC_CPR *)pc->data;
  PetscBool ascii;

  PetscFunctionBegin;
  PetscCall(PetscObjectTypeCompare((PetscObject)viewer, PETSCVIEWERASCII, &ascii));
  if (!ascii) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(PetscViewerASCIIPrintf(viewer, "  CPR: block size %" PetscInt_FMT ", pressure field %" PetscInt_FMT ", weights %s, %s\n", cpr->lbs ? cpr->lbs : cpr->bs, cpr->pfield,
                                   PCCPRWeightTypes[cpr->wtype], PCCPRCompositeTypes[cpr->ctype]));
  if (cpr->nfallback) PetscCall(PetscViewerASCIIPrintf(viewer, "  CPR: %" PetscInt_FMT " nodes fell back to raw weights\n", cpr->nfallback));
  if (cpr->kspP) {
    PetscCall(PetscViewerASCIIPrintf(viewer, "  CPR pressure stage (prefix -pc_cpr_pressure_):\n"));
    PetscCall(PetscViewerASCIIPushTab(viewer));
    PetscCall(KSPView(cpr->kspP, viewer));
    PetscCall(PetscViewerASCIIPopTab(viewer));
  }
  if (cpr->kspG) {
    PetscCall(PetscViewerASCIIPrintf(viewer, "  CPR global stage (prefix -pc_cpr_global_):\n"));
    PetscCall(PetscViewerASCIIPushTab(viewer));
    PetscCall(KSPView(cpr->kspG, viewer));
    PetscCall(PetscViewerASCIIPopTab(viewer));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PCCPRSetBlockSize_CPR(PC pc, PetscInt bs)
{
  PetscFunctionBegin;
  ((PC_CPR *)pc->data)->bs = bs;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PCCPRSetPressureField_CPR(PC pc, PetscInt f)
{
  PetscFunctionBegin;
  ((PC_CPR *)pc->data)->pfield = f;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PCCPRGetPressureKSP_CPR(PC pc, KSP *ksp)
{
  PetscFunctionBegin;
  PetscCall(CPRCreateKSPs(pc));
  *ksp = ((PC_CPR *)pc->data)->kspP;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PCCPRGetGlobalKSP_CPR(PC pc, KSP *ksp)
{
  PetscFunctionBegin;
  PetscCall(CPRCreateKSPs(pc));
  *ksp = ((PC_CPR *)pc->data)->kspG;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PCDestroy_CPR(PC pc)
{
  PC_CPR *cpr = (PC_CPR *)pc->data;

  PetscFunctionBegin;
  PetscCall(PCReset_CPR(pc));
  PetscCall(KSPDestroy(&cpr->kspP));
  PetscCall(KSPDestroy(&cpr->kspG));
  PetscCall(PetscObjectComposeFunction((PetscObject)pc, "PCCPRSetBlockSize_C", NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pc, "PCCPRSetPressureField_C", NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pc, "PCCPRGetPressureKSP_C", NULL));
  PetscCall(PetscObjectComposeFunction((PetscObject)pc, "PCCPRGetGlobalKSP_C", NULL));
  PetscCall(PetscFree(pc->data));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode PCCreate_CPR(PC pc)
{
  PC_CPR *cpr;

  PetscFunctionBegin;
  PetscCall(PetscNew(&cpr));
  cpr->pfield = 0;
  cpr->wtype  = PC_CPR_WEIGHT_QIMPES;
  cpr->ctype  = PC_CPR_MULTIPLICATIVE;
  pc->data    = (void *)cpr;

  pc->ops->setup          = PCSetUp_CPR;
  pc->ops->apply          = PCApply_CPR;
  pc->ops->reset          = PCReset_CPR;
  pc->ops->destroy        = PCDestroy_CPR;
  pc->ops->setfromoptions = PCSetFromOptions_CPR;
  pc->ops->view           = PCView_CPR;

  PetscCall(PetscObjectComposeFunction((PetscObject)pc, "PCCPRSetBlockSize_C", PCCPRSetBlockSize_CPR));
  PetscCall(PetscObjectComposeFunction((PetscObject)pc, "PCCPRSetPressureField_C", PCCPRSetPressureField_CPR));
  PetscCall(PetscObjectComposeFunction((PetscObject)pc, "PCCPRGetPressureKSP_C", PCCPRGetPressureKSP_CPR));
  PetscCall(PetscObjectComposeFunction((PetscObject)pc, "PCCPRGetGlobalKSP_C", PCCPRGetGlobalKSP_CPR));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* --------------------------------------------------------------- public API */

PetscErrorCode PCCPRSetBlockSize(PC pc, PetscInt bs)
{
  PetscFunctionBegin;
  PetscTryMethod(pc, "PCCPRSetBlockSize_C", (PC, PetscInt), (pc, bs));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PCCPRSetPressureField(PC pc, PetscInt f)
{
  PetscFunctionBegin;
  PetscTryMethod(pc, "PCCPRSetPressureField_C", (PC, PetscInt), (pc, f));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PCCPRGetPressureKSP(PC pc, KSP *ksp)
{
  PetscFunctionBegin;
  *ksp = NULL;
  PetscTryMethod(pc, "PCCPRGetPressureKSP_C", (PC, KSP *), (pc, ksp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PCCPRGetGlobalKSP(PC pc, KSP *ksp)
{
  PetscFunctionBegin;
  *ksp = NULL;
  PetscTryMethod(pc, "PCCPRGetGlobalKSP_C", (PC, KSP *), (pc, ksp));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PCCPRSetNodeDofs(PC pc, PetscInt bs, PetscInt nnodes, const PetscInt dofs[])
{
  IS is;

  PetscFunctionBegin;
  PetscCall(ISCreateGeneral(PetscObjectComm((PetscObject)pc), bs * nnodes, dofs, PETSC_COPY_VALUES, &is));
  PetscCall(ISSetBlockSize(is, bs));
  PetscCall(PetscObjectCompose((PetscObject)pc, "PCCPR_NodeDofs", (PetscObject)is));
  PetscCall(ISDestroy(&is));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PCCPRRegister(void)
{
  static PetscBool done = PETSC_FALSE;

  PetscFunctionBegin;
  if (done) PetscFunctionReturn(PETSC_SUCCESS);
  done = PETSC_TRUE;
  PetscCall(PCRegister("cpr", PCCreate_CPR));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* entry point used by PETSc when loaded with -dll_append libpccpr.so */
PETSC_EXTERN PetscErrorCode PetscDLLibraryRegister_pccpr(void)
{
  return PCCPRRegister();
}
