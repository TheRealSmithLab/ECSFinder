#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "utils.h"
#include "fold_vars.h"
#include "energy_par.h"
#include "constraints.h"
#include "exterior_loops.h"
#include "gquad.h"
#include "structured_domains.h"
#include "unstructured_domains.h"
#include "multibranch_loops.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

struct default_data {
  int                       *idx;
  unsigned char             *mx;
  unsigned char             **mx_window;
  int                       cp;
  int                       *hc_up;
  void                      *hc_dat;
  vrna_callback_hc_evaluate *hc_f;
};


/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE int
E_mb_loop_fast(vrna_fold_compound_t *vc,
               int                  i,
               int                  j,
               int                  *dmli1,
               int                  *dmli2);


PRIVATE int
E_mb_loop_fast_window(vrna_fold_compound_t  *vc,
                      int                   i,
                      int                   j,
                      int                   *dmli1,
                      int                   *dmli2);


PRIVATE int
E_mb_loop_fast_comparative(vrna_fold_compound_t *vc,
                           int                  i,
                           int                  j,
                           int                  *dmli1,
                           int                  *dmli2);


PRIVATE int
E_mb_loop_fast_comparative_window(vrna_fold_compound_t  *vc,
                                  int                   i,
                                  int                   j,
                                  int                   *dmli1,
                                  int                   *dmli2);


PRIVATE FLT_OR_DBL
exp_E_mb_loop_fast(vrna_fold_compound_t *vc,
                   int                  i,
                   int                  j,
                   FLT_OR_DBL           *qqm1);


PRIVATE FLT_OR_DBL
exp_E_mb_loop_fast_window(vrna_fold_compound_t  *vc,
                          int                   i,
                          int                   j,
                          FLT_OR_DBL            *qqm1);


PRIVATE FLT_OR_DBL
exp_E_mb_loop_fast_comparative(vrna_fold_compound_t *vc,
                               int                  i,
                               int                  j,
                               FLT_OR_DBL           *qqm1);


PRIVATE int
E_ml_stems_fast(vrna_fold_compound_t  *vc,
                int                   i,
                int                   j,
                int                   *fmi,
                int                   *dmli);


PRIVATE int
E_ml_stems_fast_comparative(vrna_fold_compound_t  *vc,
                            int                   i,
                            int                   j,
                            int                   *fmi,
                            int                   *dmli);


PRIVATE int
E_ml_stems_fast_comparative_window(vrna_fold_compound_t *vc,
                                   int                  i,
                                   int                  j,
                                   int                  *fmi,
                                   int                  *dmli);


PRIVATE int
E_ml_stems_fast_window(vrna_fold_compound_t *vc,
                       int                  i,
                       int                  j,
                       int                  *fmi,
                       int                  *dmli);


PRIVATE int
extend_fm_3p(int                  i,
             int                  j,
             int                  *fm,
             vrna_fold_compound_t *vc);


PRIVATE int
extend_fm_3p_window(int                   i,
                    int                   j,
                    int                   **fm,
                    vrna_fold_compound_t  *vc);


PRIVATE int E_mb_loop_stack(vrna_fold_compound_t  *vc,
                            int                   i,
                            int                   j);


PRIVATE int E_mb_loop_stack_window(vrna_fold_compound_t *vc,
                                   int                  i,
                                   int                  j);


PRIVATE int
BT_mb_loop(vrna_fold_compound_t *vc,
           int                  *i,
           int                  *j,
           int                  *k,
           int                  en,
           int                  *component1,
           int                  *component2);


PRIVATE int
BT_mb_loop_comparative(vrna_fold_compound_t *vc,
                       int                  *i,
                       int                  *j,
                       int                  *k,
                       int                  en,
                       int                  *component1,
                       int                  *component2);


PRIVATE int
BT_mb_loop_window(vrna_fold_compound_t  *vc,
                  int                   *i,
                  int                   *j,
                  int                   *k,
                  int                   en,
                  int                   *component1,
                  int                   *component2);


PRIVATE int
BT_mb_loop_window_comparative(vrna_fold_compound_t  *vc,
                              int                   *i,
                              int                   *j,
                              int                   *k,
                              int                   en,
                              int                   *component1,
                              int                   *component2);


PRIVATE int
BT_mb_loop_split(vrna_fold_compound_t *vc,
                 int                  *i,
                 int                  *j,
                 int                  *k,
                 int                  *l,
                 int                  *component1,
                 int                  *component2,
                 vrna_bp_stack_t      *bp_stack,
                 int                  *stack_count);


PRIVATE int
BT_mb_loop_split_comparative(vrna_fold_compound_t *vc,
                             int                  *i,
                             int                  *j,
                             int                  *k,
                             int                  *l,
                             int                  *component1,
                             int                  *component2,
                             vrna_bp_stack_t      *bp_stack,
                             int                  *stack_count);


PRIVATE int
BT_mb_loop_split_window(vrna_fold_compound_t  *vc,
                        int                   *i,
                        int                   *j,
                        int                   *k,
                        int                   *l,
                        int                   *component1,
                        int                   *component2,
                        vrna_bp_stack_t       *bp_stack,
                        int                   *stack_count);


PRIVATE int
BT_mb_loop_split_window_comparative(vrna_fold_compound_t  *vc,
                                    int                   *i,
                                    int                   *j,
                                    int                   *k,
                                    int                   *l,
                                    int                   *component1,
                                    int                   *component2,
                                    vrna_bp_stack_t       *bp_stack,
                                    int                   *stack_count);


PRIVATE unsigned char
hc_default(int            i,
           int            j,
           int            k,
           int            l,
           unsigned char  d,
           void           *data);


PRIVATE unsigned char
hc_default_window(int           i,
                  int           j,
                  int           k,
                  int           l,
                  unsigned char d,
                  void          *data);


PRIVATE unsigned char
hc_default_ext(int            i,
               int            j,
               int            k,
               int            l,
               unsigned char  d,
               void           *data);


PRIVATE unsigned char
hc_default_user(int           i,
                int           j,
                int           k,
                int           l,
                unsigned char d,
                void          *data);


PRIVATE unsigned char
hc_default_user_window(int            i,
                       int            j,
                       int            k,
                       int            l,
                       unsigned char  d,
                       void           *data);


PRIVATE unsigned char
hc_default_user_ext(int           i,
                    int           j,
                    int           k,
                    int           l,
                    unsigned char d,
                    void          *data);


PRIVATE FLT_OR_DBL
exp_E_ml_fast(vrna_fold_compound_t  *vc,
              int                   i,
              int                   j,
              vrna_mx_pf_aux_ml_t   *aux_mx);


PRIVATE FLT_OR_DBL
exp_E_ml_fast_window(vrna_fold_compound_t *vc,
                     int                  i,
                     int                  j,
                     vrna_mx_pf_aux_ml_t  *aux_mx);


PRIVATE FLT_OR_DBL
exp_E_ml_fast_comparative(vrna_fold_compound_t  *vc,
                          int                   i,
                          int                   j,
                          vrna_mx_pf_aux_ml_t   *aux_mx);


PRIVATE INLINE int get_pair_type_md(int       i,
                                    int       j,
                                    vrna_md_t *md);


PRIVATE INLINE int
get_pair_type(int   ij,
              char  *ptype);


PRIVATE INLINE int
get_pair_type_window(int  i,
                     int  j,
                     char **ptype);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC int
vrna_E_mb_loop_fast(vrna_fold_compound_t  *vc,
                    int                   i,
                    int                   j,
                    int                   *dmli1,
                    int                   *dmli2)
{
  int e = INF;

  if (vc) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          e = E_mb_loop_fast_window(vc, i, j, dmli1, dmli2);
        else
          e = E_mb_loop_fast(vc, i, j, dmli1, dmli2);

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          e = E_mb_loop_fast_comparative_window(vc, i, j, dmli1, dmli2);
        else
          e = E_mb_loop_fast_comparative(vc, i, j, dmli1, dmli2);

        break;
    }
  }

  return e;
}


PUBLIC int
vrna_E_ml_stems_fast(vrna_fold_compound_t *vc,
                     int                  i,
                     int                  j,
                     int                  *fmi,
                     int                  *dmli)
{
  int e = INF;

  if (vc) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          e = E_ml_stems_fast_window(vc, i, j, fmi, dmli);
        else
          e = E_ml_stems_fast(vc, i, j, fmi, dmli);

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          e = E_ml_stems_fast_comparative_window(vc, i, j, fmi, dmli);
        else
          e = E_ml_stems_fast_comparative(vc, i, j, fmi, dmli);

        break;
    }
  }

  return e;
}


PUBLIC FLT_OR_DBL
vrna_exp_E_mb_loop_fast(vrna_fold_compound_t  *vc,
                        int                   i,
                        int                   j,
                        FLT_OR_DBL            *qqm1)
{
  FLT_OR_DBL q = 0.;

  if (vc) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          q = exp_E_mb_loop_fast_window(vc, i, j, qqm1);
        else
          q = exp_E_mb_loop_fast(vc, i, j, qqm1);

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        q = exp_E_mb_loop_fast_comparative(vc, i, j, qqm1);
        break;
    }
  }

  return q;
}


PUBLIC int
vrna_E_mb_loop_stack(vrna_fold_compound_t *vc,
                     int                  i,
                     int                  j)
{
  int e = INF;

  if (vc) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          e = E_mb_loop_stack_window(vc, i, j);
        else
          e = E_mb_loop_stack(vc, i, j);

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        break;
    }
  }

  return e;
}


PUBLIC int
vrna_BT_mb_loop_split(vrna_fold_compound_t  *vc,
                      int                   *i,
                      int                   *j,
                      int                   *k,
                      int                   *l,
                      int                   *c1,
                      int                   *c2,
                      vrna_bp_stack_t       *bp_stack,
                      int                   *stack_count)
{
  if (vc) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          return BT_mb_loop_split_window(vc, i, j, k, l, c1, c2, bp_stack, stack_count);
        else
          return BT_mb_loop_split(vc, i, j, k, l, c1, c2, bp_stack, stack_count);

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          return BT_mb_loop_split_window_comparative(vc, i, j, k, l, c1, c2, bp_stack, stack_count);

        else
          return BT_mb_loop_split_comparative(vc, i, j, k, l, c1, c2, bp_stack, stack_count);

        break;
    }
  }

  return 0;
}


PUBLIC int
vrna_BT_mb_loop(vrna_fold_compound_t  *vc,
                int                   *i,
                int                   *j,
                int                   *k,
                int                   en,
                int                   *c1,
                int                   *c2)
{
  if (vc) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          return BT_mb_loop_window(vc, i, j, k, en, c1, c2);
        else
          return BT_mb_loop(vc, i, j, k, en, c1, c2);

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          return BT_mb_loop_window_comparative(vc, i, j, k, en, c1, c2);
        else
          return BT_mb_loop_comparative(vc, i, j, k, en, c1, c2);

        break;
    }
  }

  return 0;
}


PRIVATE INLINE int
get_pair_type_md(int        i,
                 int        j,
                 vrna_md_t  *md)
{
  int tt = md->pair[i][j];

  return (tt == 0) ? 7 : tt;
}


PRIVATE INLINE int
get_pair_type(int   ij,
              char  *ptype)
{
  int tt = (int)ptype[ij];

  return (tt == 0) ? 7 : tt;
}


PRIVATE INLINE int
get_pair_type_window(int  i,
                     int  j,
                     char **ptype)
{
  int tt = (int)ptype[i][j - i];

  return (tt == 0) ? 7 : tt;
}


PRIVATE int
E_mb_loop_fast_comparative(vrna_fold_compound_t *vc,
                           int                  i,
                           int                  j,
                           int                  *dmli1,
                           int                  *dmli2)
{
  short                     **S, **S5, **S3;
  int                       *indx, e, decomp, s, n_seq, dangle_model, tt;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n_seq         = vc->n_seq;
  indx          = vc->jindx;
  P             = vc->params;
  md            = &(P->model_details);
  hc            = vc->hc;
  scs           = vc->scs;
  dangle_model  = md->dangles;
  e             = INF;

  hc_dat_local.idx        = vc->jindx;
  hc_dat_local.mx         = hc->matrix;
  hc_dat_local.mx_window  = hc->matrix_local;
  hc_dat_local.hc_up      = hc->up_ml;
  hc_dat_local.cp         = vc->cutpoint;

  if (hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default;
  }

  /* multi-loop decomposition ------------------------*/
  if (evaluate(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
    decomp = dmli1[j - 1];

    S   = vc->S;
    S5  = vc->S5;       /* S5[s][i] holds next base 5' of i in sequence s */
    S3  = vc->S3;       /* Sl[s][i] holds next base 3' of i in sequence s */

    if (dangle_model) {
      for (s = 0; s < n_seq; s++) {
        tt      = get_pair_type_md(S[s][j], S[s][i], md);
        decomp  += E_MLstem(tt, S5[s][j], S3[s][i], P);
      }
    } else {
      for (s = 0; s < n_seq; s++) {
        tt      = get_pair_type_md(S[s][j], S[s][i], md);
        decomp  += E_MLstem(tt, -1, -1, P);
      }
    }

    if (scs) {
      for (s = 0; s < n_seq; s++) {
        if (scs[s])
          if (scs[s]->energy_bp)
            decomp += scs[s]->energy_bp[indx[j] + i];
      }
    }

    e = decomp + n_seq * P->MLclosing;
  }

  return e;
}


PRIVATE int
E_mb_loop_fast_comparative_window(vrna_fold_compound_t  *vc,
                                  int                   i,
                                  int                   j,
                                  int                   *dmli1,
                                  int                   *dmli2)
{
  short                     **S, **S5, **S3;
  int                       e, decomp, s, n_seq, dangle_model, tt;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n_seq         = vc->n_seq;
  P             = vc->params;
  md            = &(P->model_details);
  hc            = vc->hc;
  scs           = vc->scs;
  dangle_model  = md->dangles;
  e             = INF;

  hc_dat_local.idx        = vc->jindx;
  hc_dat_local.mx         = hc->matrix;
  hc_dat_local.mx_window  = hc->matrix_local;
  hc_dat_local.hc_up      = hc->up_ml;
  hc_dat_local.cp         = vc->cutpoint;

  if (hc->f) {
    evaluate            = &hc_default_user_window;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default_window;
  }

  /* multi-loop decomposition ------------------------*/
  if (evaluate(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
    decomp = dmli1[j - 1 - (i + 1)];

    S   = vc->S;
    S5  = vc->S5;       /* S5[s][i] holds next base 5' of i in sequence s */
    S3  = vc->S3;       /* Sl[s][i] holds next base 3' of i in sequence s */

    if (dangle_model) {
      for (s = 0; s < n_seq; s++) {
        tt      = get_pair_type_md(S[s][j], S[s][i], md);
        decomp  += E_MLstem(tt, S5[s][j], S3[s][i], P);
      }
    } else {
      for (s = 0; s < n_seq; s++) {
        tt      = get_pair_type_md(S[s][j], S[s][i], md);
        decomp  += E_MLstem(tt, -1, -1, P);
      }
    }

    if (scs) {
      for (s = 0; s < n_seq; s++) {
        if (scs[s])
          if (scs[s]->f)
            decomp += scs[s]->f(i, j, i + 2, j - 1, VRNA_DECOMP_PAIR_ML, scs[s]->data);
      }
    }

    e = decomp + n_seq * P->MLclosing;
  }

  return e;
}


PRIVATE int
E_mb_loop_fast(vrna_fold_compound_t *vc,
               int                  i,
               int                  j,
               int                  *dmli1,
               int                  *dmli2)
{
  short                     S_i1, S_j1, *S, *S2;
  unsigned int              *sn;
  int                       decomp, en, e, cp, *indx, *fc, ij, dangle_model, tt;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  cp            = vc->cutpoint;
  S             = vc->sequence_encoding;
  S2            = vc->sequence_encoding2;
  indx          = vc->jindx;
  sn            = vc->strand_number;
  hc            = vc->hc;
  sc            = vc->sc;
  fc            = vc->matrices->fc;
  P             = vc->params;
  md            = &(P->model_details);
  dangle_model  = md->dangles;

  /* init values */
  e       = INF;
  decomp  = INF;

  hc_dat_local.idx        = vc->jindx;
  hc_dat_local.mx         = hc->matrix;
  hc_dat_local.mx_window  = hc->matrix_local;
  hc_dat_local.hc_up      = hc->up_ml;
  hc_dat_local.cp         = vc->cutpoint;

  ij  = indx[j] + i;
  tt  = get_pair_type_md(S2[j], S2[i], md);

  if (hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default;
  }

  if (cp < 0) {
    S_i1  = S[i + 1];
    S_j1  = S[j - 1];
  } else {
    S_i1  = (sn[i] == sn[i + 1]) ? S[i + 1] : -1;
    S_j1  = (sn[j - 1] == sn[j]) ? S[j - 1] : -1;
  }

  if ((S_i1 >= 0) && (S_j1 >= 0)) {
    /* regular multi branch loop */
    /* new closing pair (i,j) with mb part [i+1,j-1] */
    if (evaluate(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
      decomp = dmli1[j - 1];

      if (decomp != INF) {
        switch (dangle_model) {
          /* no dangles */
          case 0:
            decomp += E_MLstem(tt, -1, -1, P);
            if (sc) {
              if (sc->energy_bp)
                decomp += sc->energy_bp[ij];

              if (sc->f)
                decomp += sc->f(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, sc->data);
            }

            break;

          /* double dangles */
          case 2:
            decomp += E_MLstem(tt, S_j1, S_i1, P);
            if (sc) {
              if (sc->energy_bp)
                decomp += sc->energy_bp[ij];

              if (sc->f)
                decomp += sc->f(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, sc->data);
            }

            break;

          /* normal dangles, aka dangles = 1 || 3 */
          default:
            decomp += E_MLstem(tt, -1, -1, P);
            if (sc) {
              if (sc->energy_bp)
                decomp += sc->energy_bp[ij];

              if (sc->f)
                decomp += sc->f(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, sc->data);
            }

            break;
        }
      }
    }

    if (dangle_model % 2) {
      /* dangles == 1 || dangles == 3 */
      /* new closing pair (i,j) with mb part [i+2,j-1] */
      if (evaluate(i, j, i + 2, j - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
        if (dmli2[j - 1] != INF) {
          en = dmli2[j - 1] +
               E_MLstem(tt, -1, S_i1, P) +
               P->MLbase;

          if (sc) {
            if (sc->energy_up)
              en += sc->energy_up[i + 1][1];

            if (sc->energy_bp)
              en += sc->energy_bp[ij];

            if (sc->f)
              en += sc->f(i, j, i + 2, j - 1, VRNA_DECOMP_PAIR_ML, sc->data);
          }

          decomp = MIN2(decomp, en);
        }
      }

      /* new closing pair (i,j) with mb part [i+2.j-2] */
      if (evaluate(i, j, i + 2, j - 2, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
        if (dmli2[j - 2] != INF) {
          en = dmli2[j - 2] +
               E_MLstem(tt, S_j1, S_i1, P) +
               2 * P->MLbase;

          if (sc) {
            if (sc->energy_up)
              en += sc->energy_up[i + 1][1] +
                    sc->energy_up[j - 1][1];

            if (sc->energy_bp)
              en += sc->energy_bp[ij];

            if (sc->f)
              en += sc->f(i, j, i + 2, j - 2, VRNA_DECOMP_PAIR_ML, sc->data);
          }

          decomp = MIN2(decomp, en);
        }
      }

      /* new closing pair (i,j) with mb part [i+1, j-2] */
      if (evaluate(i, j, i + 1, j - 2, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
        if (dmli1[j - 2] != INF) {
          en = dmli1[j - 2] +
               E_MLstem(tt, S_j1, -1, P) +
               P->MLbase;
          if (sc) {
            if (sc->energy_up)
              en += sc->energy_up[j - 1][1];

            if (sc->energy_bp)
              en += sc->energy_bp[ij];

            if (sc->f)
              en += sc->f(i, j, i + 1, j - 2, VRNA_DECOMP_PAIR_ML, sc->data);
          }

          decomp = MIN2(decomp, en);
        }
      }
    } /* end if dangles % 2 */

    if (decomp != INF)
      e = decomp + P->MLclosing;
  } /* end regular multibranch loop */

  if (sn[i] != sn[j]) {
    /* multibrach like cofold structure with cut somewhere between i and j */
    if (evaluate(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
      if ((fc[i + 1] != INF) && (fc[j - 1] != INF)) {
        decomp = fc[i + 1] +
                 fc[j - 1];
        switch (dangle_model) {
          case 0:
            decomp += E_ExtLoop(tt, -1, -1, P);
            break;

          case 2:
            decomp += E_ExtLoop(tt, S_j1, S_i1, P);
            break;

          default:
            decomp += E_ExtLoop(tt, -1, -1, P);
            break;
        }
      }
    }

    if (dangle_model % 2) {
      /* dangles == 1 || dangles == 3 */
      if (evaluate(i, j, i + 2, j - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
        if ((fc[i + 2] != INF) && (fc[j - 1] != INF)) {
          en = fc[i + 2] +
               fc[j - 1] +
               E_ExtLoop(tt, -1, S_i1, P);
          decomp = MIN2(decomp, en);
        }
      }

      if (evaluate(i, j, i + 1, j - 2, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
        if ((fc[i + 1] != INF) && (fc[j - 2] != INF)) {
          en = fc[i + 1] +
               fc[j - 2] +
               E_ExtLoop(tt, S_j1, -1, P);
          decomp = MIN2(decomp, en);
        }
      }

      if (evaluate(i, j, i + 2, j - 2, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
        if ((fc[i + 2] != INF) && (fc[j - 2] != INF)) {
          en = fc[i + 2] +
               fc[j - 2] +
               E_ExtLoop(tt, S_j1, S_i1, P);
          decomp = MIN2(decomp, en);
        }
      }
    }

    e = MIN2(e, decomp);
  }

  return e;
}


PRIVATE int
E_mb_loop_fast_window(vrna_fold_compound_t  *vc,
                      int                   i,
                      int                   j,
                      int                   *dmli1,
                      int                   *dmli2)
{
  short                     S_i1, S_j1, *S, *S2;
  int                       decomp, en, e, dangle_model, tt;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  S             = vc->sequence_encoding;
  S2            = vc->sequence_encoding2;
  hc            = vc->hc;
  sc            = vc->sc;
  P             = vc->params;
  md            = &(P->model_details);
  dangle_model  = md->dangles;

  /* init values */
  e       = INF;
  decomp  = INF;

  hc_dat_local.idx        = vc->jindx;
  hc_dat_local.mx         = hc->matrix;
  hc_dat_local.mx_window  = hc->matrix_local;
  hc_dat_local.hc_up      = hc->up_ml;
  hc_dat_local.cp         = vc->cutpoint;

  if (hc->f) {
    evaluate            = &hc_default_user_window;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default_window;
  }

  S_i1  = S[i + 1];
  S_j1  = S[j - 1];

  tt = get_pair_type_md(S2[j], S2[i], md);

  /* new closing pair (i,j) with mb part [i+1,j-1] */
  if (evaluate(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
    decomp = dmli1[j - 1 - (i + 1)];

    if (decomp < INF) {
      switch (dangle_model) {
        /* no dangles */
        case 0:
          decomp += E_MLstem(tt, -1, -1, P);
          break;
        /* double dangles */
        case 2:
          decomp += E_MLstem(tt, S_j1, S_i1, P);
          break;
        /* normal dangles, aka dangle_model = 1 */
        default:
          decomp += E_MLstem(tt, -1, -1, P);
          break;
      }
      if (sc) {
        if (sc->energy_bp_local)
          decomp += sc->energy_bp_local[i][j - i];

        if (sc->f)
          decomp += sc->f(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, sc->data);
      }
    }
  }

  if (dangle_model % 2) {
    /* dangles == 1 || dangles == 3 */
    /* new closing pair (i,j) with mb part [i+2,j-1] */
    if (evaluate(i, j, i + 2, j - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
      if (dmli2[j - 1 - (i + 2)] != INF) {
        en = dmli2[j - 1 - (i + 2)] +
             E_MLstem(tt, -1, S_i1, P) +
             P->MLbase;

        if (sc) {
          if (sc->energy_up)
            decomp += sc->energy_up[i + 1][1];

          if (sc->energy_bp_local)
            decomp += sc->energy_bp_local[i][j - i];

          if (sc->f)
            en += sc->f(i, j, i + 2, j - 1, VRNA_DECOMP_PAIR_ML, sc->data);
        }

        decomp = MIN2(decomp, en);
      }
    }

    /* new closing pair (i,j) with mb part [i+2.j-2] */
    if (evaluate(i, j, i + 2, j - 2, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
      if (dmli2[j - 2 - (i + 2)] != INF) {
        en = dmli2[j - 2 - (i + 2)] +
             E_MLstem(tt, S_j1, S_i1, P) +
             2 * P->MLbase;

        if (sc) {
          if (sc->energy_up)
            decomp += sc->energy_up[i + 1][1] +
                      sc->energy_up[j - 1][1];

          if (sc->energy_bp_local)
            decomp += sc->energy_bp_local[i][j - i];

          if (sc->f)
            en += sc->f(i, j, i + 2, j - 2, VRNA_DECOMP_PAIR_ML, sc->data);
        }

        decomp = MIN2(decomp, en);
      }
    }

    /* new closing pair (i,j) with mb part [i+1, j-2] */
    if (evaluate(i, j, i + 1, j - 2, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
      if (dmli1[j - 2 - (i + 1)] != INF) {
        en = dmli1[j - 2 - (i + 1)] +
             E_MLstem(tt, S_j1, -1, P) +
             P->MLbase;

        if (sc) {
          if (sc->energy_up)
            decomp += sc->energy_up[j - 1][1];

          if (sc->energy_bp_local)
            decomp += sc->energy_bp_local[i][j - i];

          if (sc->f)
            en += sc->f(i, j, i + 1, j - 2, VRNA_DECOMP_PAIR_ML, sc->data);
        }

        decomp = MIN2(decomp, en);
      }
    }
  } /* end if dangles % 2 */

  if (decomp != INF)
    e = decomp + P->MLclosing;

  return e;
}


PRIVATE int
E_mb_loop_stack(vrna_fold_compound_t  *vc,
                int                   i,
                int                   j)
{
  char                      *ptype;
  int                       *c, *fML, e, decomp, en, i1k, k1j1, ij, k, *indx, turn,
                            type, type_2, *rtype;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  indx  = vc->jindx;
  hc    = vc->hc;
  P     = vc->params;
  md    = &(P->model_details);
  turn  = md->min_loop_size;
  rtype = &(md->rtype[0]);
  c     = vc->matrices->c;
  fML   = vc->matrices->fML;
  ptype = vc->ptype;
  ij    = indx[j] + i;
  type  = get_pair_type(ij, ptype);
  sc    = vc->sc;
  e     = INF;

  hc_dat_local.idx        = vc->jindx;
  hc_dat_local.mx         = hc->matrix;
  hc_dat_local.mx_window  = hc->matrix_local;
  hc_dat_local.hc_up      = hc->up_ml;
  hc_dat_local.cp         = vc->cutpoint;


  if (hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default;
  }

  if (evaluate(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
    decomp  = INF;
    k1j1    = indx[j - 1] + i + 2 + turn + 1;
    for (k = i + 2 + turn; k < j - 2 - turn; k++, k1j1++) {
      i1k = indx[k] + i + 1;

      if (evaluate(i, j, i + 1, k, VRNA_DECOMP_ML_COAXIAL, &hc_dat_local)) {
        type_2 = rtype[get_pair_type(i1k, ptype)];

        en = c[i1k] +
             P->stack[type][type_2] +
             fML[k1j1];

        if (sc)
          if (sc->f)
            en += sc->f(i, j, i + 1, k, VRNA_DECOMP_ML_COAXIAL, sc->data);

        decomp = MIN2(decomp, en);
      }

      if (evaluate(i, j, k + 1, j - 1, VRNA_DECOMP_ML_COAXIAL, &hc_dat_local)) {
        type_2 = rtype[get_pair_type(k1j1, ptype)];

        en = c[k1j1] +
             P->stack[type][type_2] +
             fML[i1k];

        if (sc)
          if (sc->f)
            en += sc->f(i, j, k + 1, j - 1, VRNA_DECOMP_ML_COAXIAL, sc->data);

        decomp = MIN2(decomp, en);
      }
    }
    /* no TermAU penalty if coax stack */
    decomp += 2 * P->MLintern[1] + P->MLclosing;
    if (sc) {
      if (sc->energy_bp)
        decomp += sc->energy_bp[ij];

      if (sc->f)
        decomp += sc->f(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, sc->data);
    }

    e = decomp;
  }

  return e;
}


PRIVATE int
E_mb_loop_stack_window(vrna_fold_compound_t *vc,
                       int                  i,
                       int                  j)
{
  char                      **ptype;
  int                       **c, **fML, e, decomp, en, k, turn, *rtype, type, type_2;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  hc    = vc->hc;
  P     = vc->params;
  md    = &(P->model_details);
  turn  = md->min_loop_size;
  rtype = &(md->rtype[0]);
  c     = vc->matrices->c_local;
  fML   = vc->matrices->fML_local;
  sc    = vc->sc;
  e     = INF;

  hc_dat_local.idx        = vc->jindx;
  hc_dat_local.mx         = hc->matrix;
  hc_dat_local.mx_window  = hc->matrix_local;
  hc_dat_local.hc_up      = hc->up_ml;
  hc_dat_local.cp         = vc->cutpoint;


  if (hc->f) {
    evaluate            = &hc_default_user_window;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default_window;
  }

  if (evaluate(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
    ptype = vc->ptype_local;
    type  = get_pair_type_window(i, j, ptype);

    decomp = INF;
    for (k = i + 2 + turn; k < j - 2 - turn; k++) {
      if (evaluate(i, j, i + 1, k, VRNA_DECOMP_ML_COAXIAL, &hc_dat_local)) {
        type_2 = rtype[get_pair_type_window(i + 1, k, ptype)];

        en = c[i + 1][k - i - 1] +
             P->stack[type][type_2] +
             fML[k + 1][j - 1 - k - 1];

        if (sc)
          if (sc->f)
            en += sc->f(i, j, i + 1, k, VRNA_DECOMP_ML_COAXIAL, sc->data);

        decomp = MIN2(decomp, en);
      }

      if (evaluate(i, j, k + 1, j - 1, VRNA_DECOMP_ML_COAXIAL, &hc_dat_local)) {
        type_2 = rtype[get_pair_type_window(k + 1, j - 1, ptype)];

        en = c[k + 1][j - 1 - k - 1] +
             P->stack[type][type_2] +
             fML[i + 1][k - i - 1];

        if (sc)
          if (sc->f)
            en += sc->f(i, j, k + 1, j - 1, VRNA_DECOMP_ML_COAXIAL, sc->data);

        decomp = MIN2(decomp, en);
      }
    }
    /* no TermAU penalty if coax stack */
    decomp += 2 * P->MLintern[1] + P->MLclosing;
    if (sc) {
      if (sc->energy_bp_local)
        decomp += sc->energy_bp_local[i][j - i];

      if (sc->f)
        decomp += sc->f(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, sc->data);
    }

    e = decomp;
  }

  return e;
}


PUBLIC int
E_ml_rightmost_stem(int                   i,
                    int                   j,
                    vrna_fold_compound_t  *vc)
{
  if ((vc) && (vc->matrices) && (vc->matrices->fM1))
    return extend_fm_3p(i, j, vc->matrices->fM1, vc);

  return INF;
}


/*
 * compose a multibranch loop part fm[i:j]
 * by either c[i,j]/ggg[i,j] or fm[i:j-1]
 *
 * This function can be used for fM and fM1
 */
PRIVATE int
extend_fm_3p(int                  i,
             int                  j,
             int                  *fm,
             vrna_fold_compound_t *vc)
{
  short                     *S;
  unsigned int              *sn;
  int                       en, length, *indx, *c, *ggg, ij, type,
                            dangle_model, with_gquad, e, u, k, cnt, with_ud;
  vrna_param_t              *P;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_ud_t                 *domains_up;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  P             = vc->params;
  length        = vc->length;
  S             = vc->sequence_encoding;
  indx          = vc->jindx;
  sn            = vc->strand_number;
  hc            = vc->hc;
  sc            = vc->sc;
  c             = vc->matrices->c;
  ggg           = vc->matrices->ggg;
  ij            = indx[j] + i;
  type          = get_pair_type(ij, vc->ptype);
  dangle_model  = P->model_details.dangles;
  with_gquad    = P->model_details.gquad;
  domains_up    = vc->domains_up;
  with_ud       = (domains_up && domains_up->energy_cb) ? 1 : 0;
  e             = INF;

  hc_dat_local.idx    = vc->jindx;
  hc_dat_local.mx     = hc->matrix;
  hc_dat_local.hc_up  = hc->up_ml;
  hc_dat_local.cp     = vc->cutpoint;

  if (hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default;
  }

  if (sn[i - 1] == sn[i]) {
    if (sn[j] == sn[j + 1]) {
      if (evaluate(i, j, i, j, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        e = c[ij];
        if (e != INF) {
          switch (dangle_model) {
            case 2:
              e += E_MLstem(type, (i == 1) ? S[length] : S[i - 1], S[j + 1], P);
              break;

            default:
              e += E_MLstem(type, -1, -1, P);
              break;
          }
          if (sc)
            if (sc->f)
              e += sc->f(i, j, i, j, VRNA_DECOMP_ML_STEM, sc->data);
        }
      }

      if (with_gquad) {
        if (sn[i] == sn[j]) {
          en  = ggg[ij] + E_MLstem(0, -1, -1, P);
          e   = MIN2(e, en);
        }
      }
    }

    if (sn[j - 1] == sn[j]) {
      if (evaluate(i, j, i, j - 1, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
        if (fm[indx[j - 1] + i] != INF) {
          en = fm[indx[j - 1] + i] +
               P->MLbase;

          if (sc) {
            if (sc->energy_up)
              en += sc->energy_up[j][1];

            if (sc->f)
              en += sc->f(i, j, i, j - 1, VRNA_DECOMP_ML_ML, sc->data);
          }

          e = MIN2(e, en);
        }
      }
    }

    if (with_ud) {
      for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
        u = domains_up->uniq_motif_size[cnt];
        k = j - u + 1;
        if ((k > i) && (sn[j - u] == sn[j])) {
          if (evaluate(i, j, i, k - 1, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
            if (fm[indx[k - 1] + i] != INF) {
              en = domains_up->energy_cb(vc,
                                         k,
                                         j,
                                         VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                         domains_up->data);
              if (en != INF) {
                en += fm[indx[k - 1] + i] +
                      u * P->MLbase;

                if (sc) {
                  if (sc->energy_up)
                    en += sc->energy_up[k][u];

                  if (sc->f)
                    en += sc->f(i, j, i, k - 1, VRNA_DECOMP_ML_ML, sc->data);
                }

                e = MIN2(e, en);
              }
            }
          }
        }
      }
    }
  }

  return e;
}


#ifdef VRNA_WITH_SSE_IMPLEMENTATION
/* SSE modular decomposition -------------------------------*/
#include <emmintrin.h>
#include <smmintrin.h>

//http://stackoverflow.com/questions/9877700/getting-max-value-in-a-m128i-vector-with-sse
int
horizontal_min_Vec4i(__m128i x)
{
  __m128i min1  = _mm_shuffle_epi32(x, _MM_SHUFFLE(0, 0, 3, 2));
  __m128i min2  = _mm_min_epi32(x, min1);
  __m128i min3  = _mm_shuffle_epi32(min2, _MM_SHUFFLE(0, 0, 0, 1));
  __m128i min4  = _mm_min_epi32(min2, min3);

  return _mm_cvtsi128_si32(min4);
}


PRIVATE int
modular_decomposition(const int i,
                      const int ij,
                      const int j,
                      const int turn,
                      const int *fmi,
                      const int *fm)
{
  int       k       = i + turn + 1;
  int       k1j     = ij + turn + 2; //indx[j] + i + 1; //indx[j] + i + turn + 2;
  const int stop    = j - 2 - turn;
  int       decomp  = INF;
  {
    const int end = 1 + stop - k;
    int       i;
    __m128i   inf = _mm_set1_epi32(INF);

    for (i = 0; i < end - 3; i += 4) {
      __m128i   a = _mm_loadu_si128((__m128i *)&fmi[k + i]);
      __m128i   b = _mm_loadu_si128((__m128i *)&fm[k1j + i]);
      __m128i   c = _mm_add_epi32(a, b);
      /* deactivate this part if you are sure to not use any hard constraints */
#if 1
      __m128i   mask1 = _mm_cmplt_epi32(a, inf);
      __m128i   mask2 = _mm_cmplt_epi32(b, inf);
      __m128i   res   = _mm_or_si128(_mm_and_si128(mask1, c),
                                     _mm_andnot_si128(mask1, a));

      res = _mm_or_si128(_mm_and_si128(mask2, res),
                         _mm_andnot_si128(mask2, b));
      const int en = horizontal_min_Vec4i(res);
#else
      const int en = horizontal_min_Vec4i(c);
#endif
      decomp = MIN2(decomp, en);
    }
    for (; i < end; i++) {
      if ((fmi[k + i] != INF) && (fm[k1j + i] != INF)) {
        const int en = fmi[k + i] + fm[k1j + i];
        decomp = MIN2(decomp, en);
      }
    }
  }

  return decomp;
}


/* End SSE modular decomposition -------------------------------*/
#endif


PRIVATE int
extend_fm_3p_window(int                   i,
                    int                   j,
                    int                   **fm,
                    vrna_fold_compound_t  *vc)
{
  short                     *S;
  unsigned int              *sn;
  int                       en, length, **c, **ggg, type,
                            dangle_model, with_gquad, e, u, k, cnt, with_ud;
  vrna_param_t              *P;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_ud_t                 *domains_up;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  P             = vc->params;
  length        = vc->length;
  S             = vc->sequence_encoding;
  sn            = vc->strand_number;
  hc            = vc->hc;
  sc            = vc->sc;
  c             = vc->matrices->c_local;
  ggg           = vc->matrices->ggg_local;
  dangle_model  = P->model_details.dangles;
  with_gquad    = P->model_details.gquad;
  domains_up    = vc->domains_up;
  with_ud       = (domains_up && domains_up->energy_cb) ? 1 : 0;
  e             = INF;


  type = get_pair_type_window(i, j, vc->ptype_local);

  hc_dat_local.idx        = vc->jindx;
  hc_dat_local.mx_window  = hc->matrix_local;
  hc_dat_local.hc_up      = hc->up_ml;
  hc_dat_local.cp         = vc->cutpoint;

  if (hc->f) {
    evaluate            = &hc_default_user_window;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default_window;
  }

  if (sn[i - 1] == sn[i]) {
    if (sn[j] == sn[j + 1]) {
      if (evaluate(i, j, i, j, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        e = c[i][j - i];
        if (e != INF) {
          switch (dangle_model) {
            case 2:
              e += E_MLstem(type, (i == 1) ? S[length] : S[i - 1], S[j + 1], P);
              break;

            default:
              e += E_MLstem(type, -1, -1, P);
              break;
          }
          if (sc)
            if (sc->f)
              e += sc->f(i, j, i, j, VRNA_DECOMP_ML_STEM, sc->data);
        }
      }

      if (with_gquad) {
        if (sn[i] == sn[j]) {
          en = ggg[i][j - i] +
               E_MLstem(0, -1, -1, P);
          e = MIN2(e, en);
        }
      }
    }

    if (sn[j - 1] == sn[j]) {
      if (evaluate(i, j, i, j - 1, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
        if (fm[i][j - 1 - i] != INF) {
          en = fm[i][j - 1 - i] +
               P->MLbase;

          if (sc) {
            if (sc->energy_up)
              en += sc->energy_up[j][1];

            if (sc->f)
              en += sc->f(i, j, i, j - 1, VRNA_DECOMP_ML_ML, sc->data);
          }

          e = MIN2(e, en);
        }
      }
    }

    if (with_ud) {
      for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
        u = domains_up->uniq_motif_size[cnt];
        k = j - u + 1;
        if ((k > i) && (sn[j - u] == sn[j])) {
          if (evaluate(i, j, i, k - 1, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
            if (fm[i][k - 1 - i] != INF) {
              en = domains_up->energy_cb(vc,
                                         k,
                                         j,
                                         VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                         domains_up->data);
              if (en != INF) {
                en += fm[i][k - 1 - i] +
                      u * P->MLbase;

                if (sc) {
                  if (sc->energy_up)
                    en += sc->energy_up[k][u];

                  if (sc->f)
                    en += sc->f(i, j, i, k - 1, VRNA_DECOMP_ML_ML, sc->data);
                }

                e = MIN2(e, en);
              }
            }
          }
        }
      }
    }
  }

  return e;
}


PRIVATE int
E_ml_stems_fast(vrna_fold_compound_t  *vc,
                int                   i,
                int                   j,
                int                   *fmi,
                int                   *dmli)
{
  char                      *ptype;
  short                     *S;
  unsigned int              *sn;
  int                       k, en, decomp, mm5, mm3, type_2, k1j, stop, length, *indx,
                            *c, *fm, ij, dangle_model, turn, type, *rtype, circular, cp, e, u,
                            cnt, with_ud;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_param_t              *P;
  vrna_ud_t                 *domains_up;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  length        = (int)vc->length;
  ptype         = vc->ptype;
  S             = vc->sequence_encoding;
  indx          = vc->jindx;
  sn            = vc->strand_number;
  hc            = vc->hc;
  sc            = vc->sc;
  c             = vc->matrices->c;
  fm            = vc->matrices->fML;
  P             = vc->params;
  ij            = indx[j] + i;
  dangle_model  = P->model_details.dangles;
  turn          = P->model_details.min_loop_size;
  type          = get_pair_type(ij, ptype);
  rtype         = &(P->model_details.rtype[0]);
  circular      = P->model_details.circ;
  cp            = vc->cutpoint;
  domains_up    = vc->domains_up;
  with_ud       = (domains_up && domains_up->energy_cb) ? 1 : 0;
  e             = INF;

  hc_dat_local.idx    = vc->jindx;
  hc_dat_local.mx     = hc->matrix;
  hc_dat_local.hc_up  = hc->up_ml;
  hc_dat_local.cp     = vc->cutpoint;

  if (hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default;
  }

  /*
   *  extension with one unpaired nucleotide at the right (3' site)
   *  or full branch of (i,j)
   */
  e = extend_fm_3p(i, j, fm, vc);

  /*
   *  extension with one unpaired nucleotide at 5' site
   *  and all other variants which are needed for odd
   *  dangle models
   */
  if (sn[i - 1] == sn[i]) {
    if (sn[i] == sn[i + 1]) {
      if (evaluate(i, j, i + 1, j, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
        if (fm[ij + 1] != INF) {
          en = fm[ij + 1] +
               P->MLbase;

          if (sc) {
            if (sc->energy_up)
              en += sc->energy_up[i][1];

            if (sc->f)
              en += sc->f(i, j, i + 1, j, VRNA_DECOMP_ML_ML, sc->data);
          }

          e = MIN2(e, en);
        }
      }
    }

    /* extension with bound ligand on 5'site */
    if (with_ud) {
      for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
        u = domains_up->uniq_motif_size[cnt];
        k = i + u - 1;
        if ((k < j) && (sn[i] == sn[k + 1])) {
          if (evaluate(i, j, k + 1, j, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
            if (fm[ij + u] != INF) {
              en = domains_up->energy_cb(vc,
                                         i,
                                         k,
                                         VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                         domains_up->data);
              if (en != INF) {
                en += fm[ij + u] +
                      u * P->MLbase;

                if (sc) {
                  if (sc->energy_up)
                    en += sc->energy_up[i][u];

                  if (sc->f)
                    en += sc->f(i, j, k + 1, j, VRNA_DECOMP_ML_ML, sc->data);
                }

                e = MIN2(e, en);
              }
            }
          }
        }
      }
    }

    if (dangle_model % 2) {
      /* dangle_model = 1 || 3 */

      mm5 = ((i > 1) || circular) ? S[i] : -1;
      mm3 = ((j < length) || circular) ? S[j] : -1;

      if (sn[i] == sn[i + 1]) {
        if (evaluate(i, j, i + 1, j, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
          if (c[ij + 1] != INF) {
            type = get_pair_type(ij + 1, ptype);

            en = c[ij + 1] +
                 E_MLstem(type, mm5, -1, P) +
                 P->MLbase;

            if (sc) {
              if (sc->energy_up)
                en += sc->energy_up[i][1];

              if (sc->f)
                en += sc->f(i, j, i + 1, j, VRNA_DECOMP_ML_STEM, sc->data);
            }

            e = MIN2(e, en);
          }
        }
      }

      if (sn[j - 1] == sn[j]) {
        if (evaluate(i, j, i, j - 1, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
          if (c[indx[j - 1] + i] != INF) {
            type = get_pair_type(indx[j - 1] + i, ptype);

            en = c[indx[j - 1] + i] +
                 E_MLstem(type, -1, mm3, P) +
                 P->MLbase;

            if (sc) {
              if (sc->energy_up)
                en += sc->energy_up[j][1];

              if (sc->f)
                en += sc->f(i, j, i, j - 1, VRNA_DECOMP_ML_STEM, sc->data);
            }

            e = MIN2(e, en);
          }
        }
      }

      if ((sn[j - 1] == sn[j]) && (sn[i] == sn[i + 1])) {
        if (evaluate(i, j, i + 1, j - 1, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
          if (c[indx[j - 1] + i + 1] != INF) {
            type = get_pair_type(indx[j - 1] + i + 1, ptype);

            en = c[indx[j - 1] + i + 1] +
                 E_MLstem(type, mm5, mm3, P) +
                 2 * P->MLbase;

            if (sc) {
              if (sc->energy_up)
                en += sc->energy_up[j][1] + sc->energy_up[i][1];

              if (sc->f)
                en += sc->f(i, j, i + 1, j - 1, VRNA_DECOMP_ML_STEM, sc->data);
            }

            e = MIN2(e, en);
          }
        }
      }
    } /* end special cases for dangles == 1 || dangles == 3 */
  }

  /* modular decomposition -------------------------------*/
  k1j   = indx[j] + i + turn + 2;
  stop  = (cp > 0) ? (cp - 1) : (j - 2 - turn);

  /* duplicated code is faster than conditions in loop */
  if (hc->f) {
    if (sc && sc->f) {
      for (decomp = INF, k = i + 1 + turn; k <= stop; k++, k1j++) {
        if ((fmi[k] != INF) && (fm[k1j] != INF) &&
            hc->f(i, j, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
          en      = fmi[k] + fm[k1j];
          en      += sc->f(i, j, k, k + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
          decomp  = MIN2(decomp, en);
        }
      }
      k++;
      k1j++;
      for (; k <= j - 2 - turn; k++, k1j++) {
        if ((fmi[k] != INF) && (fm[k1j] != INF) &&
            hc->f(i, j, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
          en      = fmi[k] + fm[k1j];
          en      += sc->f(i, j, k, k + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
          decomp  = MIN2(decomp, en);
        }
      }
    } else {
      for (decomp = INF, k = i + 1 + turn; k <= stop; k++, k1j++) {
        if ((fmi[k] != INF) && (fm[k1j] != INF) &&
            hc->f(i, j, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
          en      = fmi[k] + fm[k1j];
          decomp  = MIN2(decomp, en);
        }
      }
      k++;
      k1j++;
      for (; k <= j - 2 - turn; k++, k1j++) {
        if ((fmi[k] != INF) && (fm[k1j] != INF) &&
            hc->f(i, j, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
          en      = fmi[k] + fm[k1j];
          decomp  = MIN2(decomp, en);
        }
      }
    }
  } else {
    if (sc && sc->f) {
      for (decomp = INF, k = i + 1 + turn; k <= stop; k++, k1j++) {
        if ((fmi[k] != INF) && (fm[k1j] != INF)) {
          en      = fmi[k] + fm[k1j];
          en      += sc->f(i, j, k, k + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
          decomp  = MIN2(decomp, en);
        }
      }
      k++;
      k1j++;
      for (; k <= j - 2 - turn; k++, k1j++) {
        if ((fmi[k] != INF) && (fm[k1j] != INF)) {
          en      = fmi[k] + fm[k1j];
          en      += sc->f(i, j, k, k + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
          decomp  = MIN2(decomp, en);
        }
      }
    } else {
#ifdef VRNA_WITH_SSE_IMPLEMENTATION

      /* modular decomposition -------------------------------*/

      decomp = modular_decomposition(i, ij, j, turn, fmi, vc->matrices->fML);
      /* end modular decomposition -------------------------------*/

#else
      for (decomp = INF, k = i + 1 + turn; k <= stop; k++, k1j++) {
        if ((fmi[k] != INF) && (fm[k1j] != INF)) {
          en      = fmi[k] + fm[k1j];
          decomp  = MIN2(decomp, en);
        }
      }
      k++;
      k1j++;
      for (; k <= j - 2 - turn; k++, k1j++) {
        if ((fmi[k] != INF) && (fm[k1j] != INF)) {
          en      = fmi[k] + fm[k1j];
          decomp  = MIN2(decomp, en);
        }
      }
#endif
    }
  }

  dmli[j] = decomp;               /* store for use in fast ML decompositon */
  e       = MIN2(e, decomp);

  /* coaxial stacking */
  if (dangle_model == 3) {
    /* additional ML decomposition as two coaxially stacked helices */
    int ik;
    k1j = indx[j] + i + turn + 2;
    for (decomp = INF, k = i + 1 + turn; k <= stop; k++, k1j++) {
      ik = indx[k] + i;
      if (evaluate(i, k, k + 1, j, VRNA_DECOMP_ML_COAXIAL_ENC, &hc_dat_local)) {
        type    = rtype[get_pair_type(ik, ptype)];
        type_2  = rtype[get_pair_type(k1j, ptype)];

        en = c[ik] +
             c[k1j] +
             P->stack[type][type_2];

        if (sc)
          if (sc->f)
            en += sc->f(i, k, k + 1, j, VRNA_DECOMP_ML_COAXIAL_ENC, sc->data);

        decomp = MIN2(decomp, en);
      }
    }
    k++;
    k1j++;
    for (; k <= j - 2 - turn; k++, k1j++) {
      ik = indx[k] + i;
      if (evaluate(i, k, k + 1, j, VRNA_DECOMP_ML_COAXIAL_ENC, &hc_dat_local)) {
        type    = rtype[get_pair_type(ik, ptype)];
        type_2  = rtype[get_pair_type(k1j, ptype)];

        en = c[ik] +
             c[k1j] +
             P->stack[type][type_2];

        if (sc)
          if (sc->f)
            en += sc->f(i, k, k + 1, j, VRNA_DECOMP_ML_COAXIAL, sc->data);

        decomp = MIN2(decomp, en);
      }
    }

    decomp += 2 * P->MLintern[1];        /* no TermAU penalty if coax stack */
#if 0
    /*
     * This is needed for Y shaped ML loops with coax stacking of
     * interior pairts, but backtracking will fail if activated
     */
    DMLi[j] = MIN2(DMLi[j], decomp);
    DMLi[j] = MIN2(DMLi[j], DMLi[j - 1] + P->MLbase);
    DMLi[j] = MIN2(DMLi[j], DMLi1[j] + P->MLbase);
    new_fML = MIN2(new_fML, DMLi[j]);
#endif
    e = MIN2(e, decomp);
  }

  fmi[j] = e;

  return e;
}


PRIVATE int
E_ml_stems_fast_window(vrna_fold_compound_t *vc,
                       int                  i,
                       int                  j,
                       int                  *fmi,
                       int                  *dmli)
{
  char                      **ptype, type, type_2, tt;
  short                     *S1;
  int                       dangle_model, **c, **fML, e, decomp, en, en2, k, turn, *rtype;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  e             = INF;
  ptype         = vc->ptype_local;
  S1            = vc->sequence_encoding;
  P             = vc->params;
  md            = &(P->model_details);
  c             = vc->matrices->c_local;
  fML           = vc->matrices->fML_local;
  hc            = vc->hc;
  sc            = vc->sc;
  type          = get_pair_type_window(i, j, ptype);
  turn          = md->min_loop_size;
  rtype         = &(md->rtype[0]);
  dangle_model  = md->dangles;

  hc_dat_local.idx        = vc->jindx;
  hc_dat_local.mx_window  = hc->matrix_local;
  hc_dat_local.hc_up      = hc->up_ml;
  hc_dat_local.cp         = vc->cutpoint;

  if (hc->f) {
    evaluate            = &hc_default_user_window;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default_window;
  }

  /*
   *  extension with one unpaired nucleotide at the right (3' site)
   *  or full branch of (i,j)
   */
  e = extend_fm_3p_window(i, j, fML, vc);

  /*
   *  extension with one unpaired nucleotide at 5' site
   *  and all other variants which are needed for odd
   *  dangle models
   */
  if (evaluate(i, j, i + 1, j, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
    decomp = fML[i + 1][j - i - 1] +
             P->MLbase;

    if (sc) {
      if (sc->energy_up)
        decomp += sc->energy_up[i][1];

      if (sc->f)
        decomp += sc->f(i, j, i + 1, j, VRNA_DECOMP_ML_ML, sc->data);
    }

    e = MIN2(e, decomp);
  }

  if (dangle_model % 2) {
    /* i+1,j */
    if (evaluate(i, j, i + 1, j, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
      tt = get_pair_type_window(i + 1, j, ptype);

      decomp = c[i + 1][j - i - 1] +
               E_MLstem(tt, S1[i], -1, P) +
               P->MLbase;

      if (sc) {
        if (sc->energy_up)
          decomp += sc->energy_up[i][1];

        if (sc->f)
          decomp += sc->f(i, j, i + 1, j, VRNA_DECOMP_ML_STEM, sc->data);
      }

      e = MIN2(e, decomp);
    }

    /* i, j-1 */
    if (evaluate(i, j, i, j - 1, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
      tt = get_pair_type_window(i, j - 1, ptype);

      decomp = c[i][j - 1 - i] +
               E_MLstem(tt, -1, S1[j], P) +
               P->MLbase;

      if (sc) {
        if (sc->energy_up)
          decomp += sc->energy_up[j][1];

        if (sc->f)
          decomp += sc->f(i, j, i, j - 1, VRNA_DECOMP_ML_STEM, sc->data);
      }

      e = MIN2(e, decomp);
    }

    /* i+1,j-1 */
    if (evaluate(i, j, i + 1, j - 1, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
      tt = get_pair_type_window(i + 1, j - 1, ptype);

      decomp = c[i + 1][j - 1 - i - 1] +
               E_MLstem(tt, S1[i], S1[j], P) +
               2 * P->MLbase;

      if (sc) {
        if (sc->energy_up)
          decomp += sc->energy_up[i][1] +
                    sc->energy_up[j][1];

        if (sc->f)
          decomp += sc->f(i, j, i + 1, j - 1, VRNA_DECOMP_ML_STEM, sc->data);
      }

      e = MIN2(e, decomp);
    }
  }

  /* modular decomposition -------------------------------*/
  if (sc && sc->f) {
    if (hc->f) {
      for (decomp = INF, k = i + 1 + turn; k <= j - 2 - turn; k++)
        if ((fmi[k - i] != INF) && (fML[k + 1][j - k - 1] != INF) &&
            hc->f(i, j, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
          en      = fmi[k - i] + fML[k + 1][j - k - 1];
          en      += sc->f(i, j, k, k + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
          decomp  = MIN2(decomp, en);
        }
    } else {
      for (decomp = INF, k = i + 1 + turn; k <= j - 2 - turn; k++) {
        en  = fmi[k - i];
        en2 = fML[k + 1][j - k - 1];
        if ((en != INF) && (en2 != INF))
          en += en2 + sc->f(i, j, k, k + 1, VRNA_DECOMP_ML_ML_ML, sc->data);

        decomp = MIN2(decomp, en);
      }
    }
  } else {
    if (hc->f) {
      for (decomp = INF, k = i + 1 + turn; k <= j - 2 - turn; k++)
        if ((fmi[k - i] != INF) && (fML[k + 1][j - k - 1] != INF) &&
            hc->f(i, j, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
          en      = fmi[k - i] + fML[k + 1][j - k - 1];
          decomp  = MIN2(decomp, en);
        }
    } else {
      for (decomp = INF, k = i + 1 + turn; k <= j - 2 - turn; k++) {
        en  = fmi[k - i];
        en2 = fML[k + 1][j - k - 1];
        if ((en != INF) && (en2 != INF))
          decomp = MIN2(decomp, en + en2);
      }
    }
  }

  dmli[j - i] = decomp;               /* store for use in ML decompositon */
  e           = MIN2(e, decomp);

  /* coaxial stacking */
  if (dangle_model == 3) {
    /* additional ML decomposition as two coaxially stacked helices */
    for (decomp = INF, k = i + 1 + turn; k <= j - 2 - turn; k++) {
      if (evaluate(i, k, k + 1, j, VRNA_DECOMP_ML_COAXIAL_ENC, &hc_dat_local)) {
        type    = rtype[get_pair_type_window(i, k, ptype)];
        type_2  = rtype[get_pair_type_window(k + 1, j, ptype)];

        en = c[i][k - i] +
             c[k + 1][j - k - 1] +
             P->stack[type][type_2];

        if (sc)
          if (sc->f)
            en += sc->f(i, k, k + 1, j, VRNA_DECOMP_ML_COAXIAL_ENC, sc->data);

        decomp = MIN2(decomp, en);
      }
    }

    decomp += 2 * P->MLintern[1];          /* no TermAU penalty if coax stack */
#if 0
    /* This is needed for Y shaped ML loops with coax stacking of
    * interior pairts, but backtracking will fail if activated */
    dmli[j - i] = MIN2(dmli[j - i], decomp);
    dmli[j - i] = MIN2(dmli[j - i], dmli[j - 1 - i] + P->MLbase);
    dmli[j - i] = MIN2(dmli[j - i], dmli1[j - (i + 1)] + P->MLbase);
    e           = MIN2(e, dmli[j - i]);
#endif
    e = MIN2(e, decomp);
  }

  fmi[j - i] = e;

  return e;
}


PRIVATE int
E_ml_stems_fast_comparative(vrna_fold_compound_t  *vc,
                            int                   i,
                            int                   j,
                            int                   *fmi,
                            int                   *dmli)
{
  unsigned char             *hard_constraints;
  short                     **S, **S5, **S3;
  unsigned int              **a2s;
  int                       e, energy, *c, *fML, *ggg, ij, *indx, s, n_seq, k,
                            dangle_model, decomp, turn, tt;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_mx_mfe_t             *matrices;
  vrna_hc_t                 *hc;
  vrna_sc_t                 **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n_seq             = vc->n_seq;
  matrices          = vc->matrices;
  P                 = vc->params;
  md                = &(P->model_details);
  c                 = matrices->c;
  fML               = matrices->fML;
  ggg               = matrices->ggg;
  indx              = vc->jindx;
  hc                = vc->hc;
  scs               = vc->scs;
  hard_constraints  = hc->matrix;
  dangle_model      = md->dangles;
  turn              = md->min_loop_size;
  a2s               = vc->a2s;
  ij                = indx[j] + i;
  e                 = INF;

  hc_dat_local.idx    = vc->jindx;
  hc_dat_local.mx     = hc->matrix;
  hc_dat_local.hc_up  = hc->up_ml;
  hc_dat_local.cp     = vc->cutpoint;

  if (hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default;
  }

  if (evaluate(i, j, i + 1, j, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
    energy = fML[ij + 1] + n_seq * P->MLbase;
    if (scs) {
      for (s = 0; s < n_seq; s++) {
        if (scs[s])
          if (scs[s]->energy_up)
            energy += scs[s]->energy_up[a2s[s][i]][1];
      }
    }

    e = MIN2(e, energy);
  }

  if (evaluate(i, j, i, j - 1, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
    energy = fML[indx[j - 1] + i] + n_seq * P->MLbase;
    if (scs) {
      for (s = 0; s < n_seq; s++) {
        if (scs[s])
          if (scs[s]->energy_up)
            energy += scs[s]->energy_up[a2s[s][j]][1];
      }
    }

    e = MIN2(e, energy);
  }

  if (evaluate(i, j, i, j, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
    energy = c[ij];

    S   = vc->S;
    S5  = vc->S5;       /* S5[s][i] holds next base 5' of i in sequence s */
    S3  = vc->S3;       /* Sl[s][i] holds next base 3' of i in sequence s */

    if (dangle_model) {
      for (s = 0; s < n_seq; s++) {
        tt      = get_pair_type_md(S[s][i], S[s][j], md);
        energy  += E_MLstem(tt, S5[s][i], S3[s][j], P);
      }
    } else {
      for (s = 0; s < n_seq; s++) {
        tt      = get_pair_type_md(S[s][i], S[s][j], md);
        energy  += E_MLstem(tt, -1, -1, P);
      }
    }

    e = MIN2(e, energy);
  }

  if (md->gquad) {
    decomp  = ggg[indx[j] + i] + n_seq * E_MLstem(0, -1, -1, P);
    e       = MIN2(e, decomp);
  }

  /* modular decomposition -------------------------------*/
  decomp = INF;
  if (hc->f) {
    for (k = i + 1 + turn; k <= j - 2 - turn; k++) {
      if (hc->f(i, j, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
        energy  = fmi[k] + fML[indx[j] + k + 1];
        decomp  = (energy < decomp) ? energy : decomp;
      }
    }
  } else {
    for (k = i + 1 + turn; k <= j - 2 - turn; k++) {
      energy  = fmi[k] + fML[indx[j] + k + 1];
      decomp  = (energy < decomp) ? energy : decomp;
    }
  }

  dmli[j] = decomp; /* store for later use in ML decompositon */

  e = MIN2(e, decomp);

  fmi[j] = e; /* store for later use in ML decompositon */

  return e;
}


PRIVATE int
E_ml_stems_fast_comparative_window(vrna_fold_compound_t *vc,
                                   int                  i,
                                   int                  j,
                                   int                  *fmi,
                                   int                  *dmli)
{
  unsigned char             **hard_constraints;
  short                     **S, **S5, **S3;
  unsigned int              **a2s;
  int                       e, energy, **c, **fML, **ggg, s, n_seq, k,
                            dangle_model, decomp, turn, tt;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_mx_mfe_t             *matrices;
  vrna_hc_t                 *hc;
  vrna_sc_t                 **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n_seq             = vc->n_seq;
  matrices          = vc->matrices;
  P                 = vc->params;
  md                = &(P->model_details);
  c                 = matrices->c_local;
  fML               = matrices->fML_local;
  ggg               = matrices->ggg_local;
  hc                = vc->hc;
  scs               = vc->scs;
  hard_constraints  = hc->matrix_local;
  dangle_model      = md->dangles;
  turn              = md->min_loop_size;
  a2s               = vc->a2s;
  e                 = INF;

  hc_dat_local.idx        = vc->jindx;
  hc_dat_local.mx_window  = hc->matrix_local;
  hc_dat_local.hc_up      = hc->up_ml;
  hc_dat_local.cp         = vc->cutpoint;

  if (hc->f) {
    evaluate            = &hc_default_user_window;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default_window;
  }

  if (evaluate(i, j, i + 1, j, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
    energy = fML[i + 1][j - (i + 1)] + n_seq * P->MLbase;
    if (scs) {
      for (s = 0; s < n_seq; s++) {
        if (scs[s])
          if (scs[s]->energy_up)
            energy += scs[s]->energy_up[a2s[s][i]][1];
      }
    }

    e = MIN2(e, energy);
  }

  if (evaluate(i, j, i, j - 1, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
    energy = fML[i][j - 1 - i] + n_seq * P->MLbase;
    if (scs) {
      for (s = 0; s < n_seq; s++) {
        if (scs[s])
          if (scs[s]->energy_up)
            energy += scs[s]->energy_up[a2s[s][j]][1];
      }
    }

    e = MIN2(e, energy);
  }

  if (evaluate(i, j, i, j, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
    energy = c[i][j - i];

    S   = vc->S;
    S5  = vc->S5;       /* S5[s][i] holds next base 5' of i in sequence s */
    S3  = vc->S3;       /* Sl[s][i] holds next base 3' of i in sequence s */

    if (dangle_model) {
      for (s = 0; s < n_seq; s++) {
        tt      = get_pair_type_md(S[s][i], S[s][j], md);
        energy  += E_MLstem(tt, S5[s][i], S3[s][j], P);
      }
    } else {
      for (s = 0; s < n_seq; s++) {
        tt      = get_pair_type_md(S[s][i], S[s][j], md);
        energy  += E_MLstem(tt, -1, -1, P);
      }
    }

    e = MIN2(e, energy);
  }

  if (md->gquad) {
    decomp = ggg[i][j - i] +
             n_seq * E_MLstem(0, -1, -1, P);
    e = MIN2(e, decomp);
  }

  /* modular decomposition -------------------------------*/
  decomp = INF;
  if (hc->f) {
    for (k = i + 1 + turn; k <= j - 2 - turn; k++) {
      if (hc->f(i, j, k, k + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
        energy  = fmi[k - i] + fML[k + 1][j - (k + 1)];
        decomp  = (decomp > energy) ? energy : decomp;
      }
    }
  } else {
    for (k = i + 1 + turn; k <= j - 2 - turn; k++) {
      energy  = fmi[k - i] + fML[k + 1][j - (k + 1)];
      decomp  = (decomp > energy) ? energy : decomp;
    }
  }

  dmli[j - i] = decomp; /* store for later use in ML decompositon */

  e = MIN2(e, decomp);

  fmi[j - i] = e; /* store for later use in ML decompositon */

  return e;
}


PRIVATE FLT_OR_DBL
exp_E_mb_loop_fast(vrna_fold_compound_t *vc,
                   int                  i,
                   int                  j,
                   FLT_OR_DBL           *qqm1)
{
  char                      *ptype;
  short                     *S1;
  unsigned int              *sn;
  int                       ij, k, kl, *my_iindx, *jindx, *rtype, tt;
  FLT_OR_DBL                qbt1, temp, qqqmmm, *qm, *scale, expMLclosing;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_exp_param_t          *pf_params;
  vrna_md_t                 *md;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  my_iindx      = vc->iindx;
  jindx         = vc->jindx;
  sc            = vc->sc;
  ptype         = vc->ptype;
  S1            = vc->sequence_encoding;
  qm            = vc->exp_matrices->qm;
  scale         = vc->exp_matrices->scale;
  pf_params     = vc->exp_params;
  md            = &(pf_params->model_details);
  ij            = jindx[j] + i;
  sn            = vc->strand_number;
  hc            = vc->hc;
  expMLclosing  = pf_params->expMLclosing;
  qbt1          = 0.;

  hc_dat_local.idx    = vc->jindx;
  hc_dat_local.mx     = hc->matrix;
  hc_dat_local.hc_up  = hc->up_ml;
  hc_dat_local.cp     = vc->cutpoint;

  if (hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default;
  }

  /* multiple stem loop contribution */
  if (evaluate(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML,
               &hc_dat_local) && (sn[i] == sn[i + 1]) && (sn[j - 1] == sn[j])) {
    rtype = &(md->rtype[0]);
    tt    = rtype[get_pair_type(ij, ptype)];

    qqqmmm = expMLclosing *
             exp_E_MLstem(tt, S1[j - 1], S1[i + 1], pf_params) *
             scale[2];

    temp  = 0.0;
    kl    = my_iindx[i + 1] - (i + 1);

    if (sc) {
      if (sc->exp_energy_bp)
        qqqmmm *= sc->exp_energy_bp[my_iindx[i] - j];

      if (sc->exp_f) {
        qqqmmm *= sc->exp_f(i, j, i, j, VRNA_DECOMP_PAIR_ML, sc->data);

        if (hc->f) {
          for (k = i + 2; k <= j - 1; k++, kl--) {
            if ((sn[k - 1] == sn[k]) &&
                (hc->f(i + 1, j - 1, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc->data))) {
              temp += qm[kl] *
                      qqm1[k] *
                      sc->exp_f(i + 1, j - 1, k - 1, k, VRNA_DECOMP_ML_ML_ML, sc->data);
            }
          }
        } else {
          for (k = i + 2; k <= j - 1; k++, kl--) {
            if (sn[k - 1] == sn[k]) {
              temp += qm[kl] *
                      qqm1[k] *
                      sc->exp_f(i + 1, j - 1, k - 1, k, VRNA_DECOMP_ML_ML_ML, sc->data);
            }
          }
        }
      } else {
        if (hc->f) {
          for (k = i + 2; k <= j - 1; k++, kl--) {
            if ((sn[k - 1] == sn[k]) &&
                (hc->f(i + 1, j - 1, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc->data)))
              temp += qm[kl] *
                      qqm1[k];
          }
        } else {
          for (k = i + 2; k <= j - 1; k++, kl--) {
            if (sn[k - 1] == sn[k])
              temp += qm[kl] *
                      qqm1[k];
          }
        }
      }
    } else {
      if (hc->f) {
        for (k = i + 2; k <= j - 1; k++, kl--) {
          if ((sn[k - 1] == sn[k]) &&
              (hc->f(i + 1, j - 1, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc->data)))
            temp += qm[kl] *
                    qqm1[k];
        }
      } else {
        for (k = i + 2; k <= j - 1; k++, kl--) {
          if (sn[k - 1] == sn[k])
            temp += qm[kl] *
                    qqm1[k];
        }
      }
    }

    qbt1 += temp * qqqmmm;
  }

  return qbt1;
}


PRIVATE FLT_OR_DBL
exp_E_mb_loop_fast_window(vrna_fold_compound_t  *vc,
                          int                   i,
                          int                   j,
                          FLT_OR_DBL            *qqm1)
{
  char                      **ptype;
  short                     *S1;
  unsigned int              *sn;
  int                       ij, k, *rtype, tt;
  FLT_OR_DBL                qbt1, temp, qqqmmm, **qm, *scale, expMLclosing;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_exp_param_t          *pf_params;
  vrna_md_t                 *md;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  sc            = vc->sc;
  ptype         = vc->ptype_local;
  S1            = vc->sequence_encoding;
  qm            = vc->exp_matrices->qm_local;
  scale         = vc->exp_matrices->scale;
  pf_params     = vc->exp_params;
  md            = &(pf_params->model_details);
  sn            = vc->strand_number;
  hc            = vc->hc;
  expMLclosing  = pf_params->expMLclosing;
  qbt1          = 0.;

  hc_dat_local.idx        = vc->jindx;
  hc_dat_local.mx_window  = hc->matrix_local;
  hc_dat_local.hc_up      = hc->up_ml;
  hc_dat_local.cp         = vc->cutpoint;

  if (hc->f) {
    evaluate            = &hc_default_user_window;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default_window;
  }

  /* multiple stem loop contribution */
  if (evaluate(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML,
               &hc_dat_local) && (sn[i] == sn[i + 1]) && (sn[j - 1] == sn[j])) {
    rtype = &(md->rtype[0]);
    tt    = rtype[get_pair_type_window(i, j + i, ptype)];

    qqqmmm = expMLclosing *
             exp_E_MLstem(tt, S1[j - 1], S1[i + 1], pf_params) *
             scale[2];

    temp = 0.0;

    if (sc) {
      if (sc->exp_energy_bp_local)
        qqqmmm *= sc->exp_energy_bp_local[i][j - i];

      if (sc->exp_f) {
        qqqmmm *= sc->exp_f(i, j, i, j, VRNA_DECOMP_PAIR_ML, sc->data);

        if (hc->f) {
          for (k = i + 2; k <= j - 1; k++) {
            if ((sn[k - 1] == sn[k]) &&
                (hc->f(i + 1, j - 1, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc->data))) {
              temp += qm[i + 1][k - 1] *
                      qqm1[k] *
                      sc->exp_f(i + 1, j - 1, k - 1, k, VRNA_DECOMP_ML_ML_ML, sc->data);
            }
          }
        } else {
          for (k = i + 2; k <= j - 1; k++) {
            if (sn[k - 1] == sn[k]) {
              temp += qm[i + 1][k - 1] *
                      qqm1[k] *
                      sc->exp_f(i + 1, j - 1, k - 1, k, VRNA_DECOMP_ML_ML_ML, sc->data);
            }
          }
        }
      } else {
        if (hc->f) {
          for (k = i + 2; k <= j - 1; k++) {
            if ((sn[k - 1] == sn[k]) &&
                (hc->f(i + 1, j - 1, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc->data)))
              temp += qm[i + 1][k - 1] *
                      qqm1[k];
          }
        } else {
          for (k = i + 2; k <= j - 1; k++) {
            if (sn[k - 1] == sn[k])
              temp += qm[i + 1][k - 1] *
                      qqm1[k];
          }
        }
      }
    } else {
      if (hc->f) {
        for (k = i + 2; k <= j - 1; k++) {
          if ((sn[k - 1] == sn[k]) &&
              (hc->f(i + 1, j - 1, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc->data)))
            temp += qm[i + 1][k - 1] *
                    qqm1[k];
        }
      } else {
        for (k = i + 2; k <= j - 1; k++) {
          if (sn[k - 1] == sn[k])
            temp += qm[i + 1][k - 1] *
                    qqm1[k];
        }
      }
    }

    qbt1 += temp * qqqmmm;
  }

  return qbt1;
}


PRIVATE FLT_OR_DBL
exp_E_mb_loop_fast_comparative(vrna_fold_compound_t *vc,
                               int                  i,
                               int                  j,
                               FLT_OR_DBL           *qqm1)
{
  short                     **S, **S5, **S3;
  int                       k, kl, *my_iindx, tt, n_seq, s;
  FLT_OR_DBL                qbt1, temp, qqqmmm, *qm, *scale, expMLclosing;
  vrna_hc_t                 *hc;
  vrna_sc_t                 **scs;
  vrna_exp_param_t          *pf_params;
  vrna_md_t                 *md;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  my_iindx      = vc->iindx;
  qm            = vc->exp_matrices->qm;
  scale         = vc->exp_matrices->scale;
  pf_params     = vc->exp_params;
  md            = &(pf_params->model_details);
  hc            = vc->hc;
  expMLclosing  = pf_params->expMLclosing;
  qbt1          = 0.;

  hc_dat_local.idx    = vc->jindx;
  hc_dat_local.mx     = hc->matrix;
  hc_dat_local.hc_up  = hc->up_ml;
  hc_dat_local.cp     = vc->cutpoint;

  if (hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default;
  }

  /* multiple stem loop contribution */
  if (evaluate(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
    S     = vc->S;
    S5    = vc->S5;       /* S5[s][i] holds next base 5' of i in sequence s */
    S3    = vc->S3;       /* Sl[s][i] holds next base 3' of i in sequence s */
    scs   = vc->scs;
    n_seq = vc->n_seq;

    qqqmmm = 1.;

    for (s = 0; s < n_seq; s++) {
      tt      = get_pair_type_md(S[s][j], S[s][i], md);
      qqqmmm  *= exp_E_MLstem(tt, S5[s][j], S3[s][i], pf_params) *
                 expMLclosing;
    }

    if (scs) {
      for (s = 0; s < n_seq; s++) {
        if (scs[s]) {
          if (scs[s]->exp_energy_bp)
            qqqmmm *= scs[s]->exp_energy_bp[my_iindx[i] - j];

          if (scs[s]->f)
            qqqmmm *= scs[s]->f(i, j, i + 1, j - 1, VRNA_DECOMP_PAIR_ML, scs[s]->data);
        }
      }
    }

    /* multi-loop loop contribution */
    temp  = 0.;
    kl    = my_iindx[i + 1] - (i + 1);

    if (hc->f) {
      if (scs) {
        for (k = i + 2; k <= j - 1; k++, kl--) {
          if (hc->f(i + 1, j - 1, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc->data)) {
            qbt1 = qm[kl] * qqm1[k];
            for (s = 0; s < n_seq; s++)
              if (scs[s] && scs[s]->f)
                qbt1 *= scs[s]->f(i + 1, j - 1, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc->data);

            temp += qbt1;
          }
        }
      } else {
        for (k = i + 2; k <= j - 1; k++, kl--)
          if (hc->f(i + 1, j - 1, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc->data))
            temp += qm[kl] * qqm1[k];
      }
    } else {
      if (scs) {
        for (k = i + 2; k <= j - 1; k++, kl--) {
          qbt1 = qm[kl] * qqm1[k];
          for (s = 0; s < n_seq; s++)
            if (scs[s] && scs[s]->f)
              qbt1 *= scs[s]->f(i + 1, j - 1, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc->data);

          temp += qbt1;
        }
      } else {
        for (k = i + 2; k <= j - 1; k++, kl--)
          temp += qm[kl] * qqm1[k];
      }
    }

    temp *= scale[2];

    qbt1 = temp * qqqmmm;
  }

  return qbt1;
}


/*
 #################################
 # Backtracking functions below  #
 #################################
 */
PUBLIC int
vrna_BT_mb_loop_fake(vrna_fold_compound_t *vc,
                     int                  *u,
                     int                  *i,
                     int                  *j,
                     vrna_bp_stack_t      *bp_stack,
                     int                  *stack_count)
{
  char                      *ptype;
  short                     mm5, mm3, *S1;
  unsigned int              *sn;
  int                       length, ii, jj, k, en, cp, fij, fi, *my_c, *my_fc, *my_ggg,
                            *idx, with_gquad, dangle_model, turn, type;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  cp            = vc->cutpoint;
  length        = vc->length;
  P             = vc->params;
  md            = &(P->model_details);
  sn            = vc->strand_number;
  hc            = vc->hc;
  sc            = vc->sc;
  S1            = vc->sequence_encoding;
  ptype         = vc->ptype;
  idx           = vc->jindx;
  my_c          = vc->matrices->c;
  my_fc         = vc->matrices->fc;
  my_ggg        = vc->matrices->ggg;
  turn          = md->min_loop_size;
  with_gquad    = md->gquad;
  dangle_model  = md->dangles;

  hc_dat_local.idx    = vc->jindx;
  hc_dat_local.mx     = hc->matrix;
  hc_dat_local.hc_up  = hc->up_ext;
  hc_dat_local.cp     = vc->cutpoint;

  if (hc->f) {
    evaluate            = &hc_default_user_ext;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default_ext;
  }

  ii  = *i;
  jj  = *j;

  if (ii < cp) {
    /* 'lower' part (fc[i<cut,j=cut-1]) */

    /* nibble off unpaired 5' bases */
    do {
      fij = my_fc[ii];
      fi  = INF;

      if (evaluate(ii, jj, ii + 1, jj, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
        fi = my_fc[ii + 1];

        if (sc)
          if (sc->energy_up)
            fi += sc->energy_up[ii][1];
      }

      if (++ii == jj)
        break;
    } while (fij == fi);
    ii--;

    if (jj < ii + turn + 2) {
      /* no more pairs */
      *u = *i = *j = -1;
      return 1;
    }

    mm5 = (ii > 1 && (sn[ii - 1] == sn[ii])) ? S1[ii - 1] : -1;

    /* i or i+1 is paired. Find pairing partner */
    switch (dangle_model) {
      case 0:
        for (k = ii + turn + 1; k <= jj; k++) {
          if (evaluate(ii, jj, k, k + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
            type = get_pair_type(idx[k] + ii, ptype);

            if (fij == my_fc[k + 1] + my_c[idx[k] + ii] + E_ExtLoop(type, -1, -1, P)) {
              bp_stack[++(*stack_count)].i  = ii;
              bp_stack[(*stack_count)].j    = k;
              *u                            = k + 1;
              *i                            = ii;
              *j                            = k;
              return 1;
            }
          }

          if (with_gquad) {
            if (fij == my_fc[k + 1] + my_ggg[idx[k] + ii]) {
              *u  = k + 1;
              *i  = *j = -1;
              vrna_BT_gquad_mfe(vc, ii, k, bp_stack, stack_count);
              return 1;
            }
          }
        }
        break;

      case 2:
        for (k = ii + turn + 1; k <= jj; k++) {
          if (evaluate(ii, jj, k, k + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
            mm3   = (sn[k] == sn[k + 1]) ? S1[k + 1] : -1;
            type  = get_pair_type(idx[k] + ii, ptype);

            if (fij == my_fc[k + 1] + my_c[idx[k] + ii] + E_ExtLoop(type, mm5, mm3, P)) {
              bp_stack[++(*stack_count)].i  = ii;
              bp_stack[(*stack_count)].j    = k;
              *u                            = k + 1;
              *i                            = ii;
              *j                            = k;
              return 1;
            }
          }

          if (with_gquad) {
            if (fij == my_fc[k + 1] + my_ggg[idx[k] + ii]) {
              *u  = k + 1;
              *i  = *j = -1;
              vrna_BT_gquad_mfe(vc, ii, k, bp_stack, stack_count);
              return 1;
            }
          }
        }
        break;

      default:
        for (k = ii + turn + 1; k <= jj; k++) {
          if (evaluate(ii, jj, k, k + 1, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
            type = get_pair_type(idx[k] + ii, ptype);

            if (fij == my_fc[k + 1] + my_c[idx[k] + ii] + E_ExtLoop(type, -1, -1, P)) {
              bp_stack[++(*stack_count)].i  = ii;
              bp_stack[(*stack_count)].j    = k;
              *u                            = k + 1;
              *i                            = ii;
              *j                            = k;
              return 1;
            }
          }

          if (evaluate(ii, jj, k, k + 2, VRNA_DECOMP_EXT_STEM_EXT, &hc_dat_local)) {
            mm3 = (sn[k] == sn[k + 1]) ? S1[k + 1] : -1;
            en  = my_c[idx[k] + ii];
            if (sc)
              if (sc->energy_up)
                en += sc->energy_up[k + 1][1];

            if (fij == my_fc[k + 2] + en + E_ExtLoop(type, -1, mm3, P)) {
              bp_stack[++(*stack_count)].i  = ii;
              bp_stack[(*stack_count)].j    = k;
              *u                            = k + 2;
              *i                            = ii;
              *j                            = k;
              return 1;
            }
          }

          if (with_gquad) {
            if (fij == my_fc[k + 1] + my_ggg[idx[k] + ii]) {
              *u  = k + 1;
              *i  = *j = -1;
              vrna_BT_gquad_mfe(vc, ii, k, bp_stack, stack_count);
              return 1;
            }
          }

          if (evaluate(ii, jj, k, k + 1, VRNA_DECOMP_EXT_STEM_EXT1, &hc_dat_local)) {
            mm5   = (sn[ii] == sn[ii + 1]) ? S1[ii] : -1;
            mm3   = (sn[k] == sn[k + 1]) ? S1[k + 1] : -1;
            type  = get_pair_type(idx[k] + ii + 1, ptype);

            en = my_c[idx[k] + ii + 1];
            if (sc)
              if (sc->energy_up)
                en += sc->energy_up[ii][1];

            if (fij == en + my_fc[k + 1] + E_ExtLoop(type, mm5, -1, P)) {
              bp_stack[++(*stack_count)].i  = ii + 1;
              bp_stack[(*stack_count)].j    = k;
              *u                            = k + 1;
              *i                            = ii + 1;
              *j                            = k;
              return 1;
            }
          }

          if ((k < jj) && (evaluate(ii, jj, k, k + 2, VRNA_DECOMP_EXT_STEM_EXT1, &hc_dat_local))) {
            mm5   = (sn[ii] == sn[ii + 1]) ? S1[ii] : -1;
            mm3   = (sn[k] == sn[k + 1]) ? S1[k + 1] : -1;
            type  = get_pair_type(idx[k] + ii + 1, ptype);

            en = my_c[idx[k] + ii + 1];
            if (sc)
              if (sc->energy_up)
                en += sc->energy_up[k + 1][1];

            if (fij == en + my_fc[k + 2] + E_ExtLoop(type, mm5, mm3, P)) {
              bp_stack[++(*stack_count)].i  = ii + 1;
              bp_stack[(*stack_count)].j    = k;
              *u                            = k + 2;
              *i                            = ii + 1;
              *j                            = k;
              return 1;
            }
          }
        }
        break;
    }
  } else {
    /* 'upper' part (fc[i=cut,j>cut]) */

    /* nibble off unpaired 3' bases */
    do {
      fij = my_fc[jj];
      fi  = INF;

      if (evaluate(ii, jj, ii, jj - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local)) {
        fi = my_fc[jj - 1];

        if (sc)
          if (sc->energy_up)
            fi += sc->energy_up[jj][1];
      }

      if (--jj == ii)
        break;
    } while (fij == fi);
    jj++;

    if (jj < ii + turn + 2) {
      /* no more pairs */
      *u = *i = *j = -1;
      return 1;
    }

    /* j or j-1 is paired. Find pairing partner */
    mm3 = ((jj < length) && (sn[jj] == sn[jj + 1])) ? S1[jj + 1] : -1;
    switch (dangle_model) {
      case 0:
        for (k = jj - turn - 1; k >= ii; k--) {
          if (with_gquad) {
            if (fij == my_fc[k - 1] + my_ggg[idx[jj] + k]) {
              *u  = k - 1;
              *i  = *j = -1;
              vrna_BT_gquad_mfe(vc, k, jj, bp_stack, stack_count);
              return 1;
            }
          }

          if (evaluate(ii, jj, k - 1, k, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
            type = get_pair_type(idx[jj] + k, ptype);

            en = my_c[idx[jj] + k];
            if (sn[k] != sn[jj])
              en += P->DuplexInit;

            if (fij == my_fc[k - 1] + en + E_ExtLoop(type, -1, -1, P)) {
              bp_stack[++(*stack_count)].i  = k;
              bp_stack[(*stack_count)].j    = jj;
              *u                            = k - 1;
              *i                            = k;
              *j                            = jj;
              return 1;
            }
          }
        }
        break;

      case 2:
        for (k = jj - turn - 1; k >= ii; k--) {
          if (with_gquad) {
            if (fij == my_fc[k - 1] + my_ggg[idx[jj] + k]) {
              *u  = k - 1;
              *i  = *j = -1;
              vrna_BT_gquad_mfe(vc, k, jj, bp_stack, stack_count);
              return 1;
            }
          }

          if (evaluate(ii, jj, k - 1, k, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
            mm5   = ((k > 1) && (sn[k - 1] == sn[k])) ? S1[k - 1] : -1;
            type  = get_pair_type(idx[jj] + k, ptype);

            en = my_c[idx[jj] + k];
            if (sn[k] != sn[jj])
              en += P->DuplexInit;

            if (fij == my_fc[k - 1] + en + E_ExtLoop(type, mm5, mm3, P)) {
              bp_stack[++(*stack_count)].i  = k;
              bp_stack[(*stack_count)].j    = jj;
              *u                            = k - 1;
              *i                            = k;
              *j                            = jj;
              return 1;
            }
          }
        }
        break;

      default:
        for (k = jj - turn - 1; k >= ii; k--) {
          if (with_gquad) {
            if (fij == my_fc[k - 1] + my_ggg[idx[jj] + k]) {
              *u  = k - 1;
              *i  = *j = -1;
              vrna_BT_gquad_mfe(vc, k, jj, bp_stack, stack_count);
              return 1;
            }
          }

          if (evaluate(ii, jj, k - 1, k, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
            type = get_pair_type(idx[jj] + k, ptype);

            en = my_c[idx[jj] + k];
            if (sn[k] != sn[jj])
              en += P->DuplexInit;

            if (fij == my_fc[k - 1] + en + E_ExtLoop(type, -1, -1, P)) {
              bp_stack[++(*stack_count)].i  = k;
              bp_stack[(*stack_count)].j    = jj;
              *u                            = k - 1;
              *i                            = k;
              *j                            = jj;
              return 1;
            }
          }

          if ((k > 1) && (sn[k - 1] == sn[k]) &&
              evaluate(ii, jj, k - 2, k, VRNA_DECOMP_EXT_EXT_STEM, &hc_dat_local)) {
            type = get_pair_type(idx[jj] + k, ptype);

            en = my_c[idx[jj] + k];
            if (sn[k] != sn[jj])
              en += P->DuplexInit;

            mm5 = S1[k - 1];

            if (sc)
              if (sc->energy_up)
                en += sc->energy_up[k - 1][1];

            if (fij == my_fc[k - 2] + en + E_ExtLoop(type, mm5, -1, P)) {
              bp_stack[++(*stack_count)].i  = k;
              bp_stack[(*stack_count)].j    = jj;
              *u                            = k - 2;
              *i                            = k;
              *j                            = jj;
              return 1;
            }
          }

          if ((sn[jj - 1] == sn[jj]) &&
              evaluate(ii, jj, k - 1, k, VRNA_DECOMP_EXT_EXT_STEM1, &hc_dat_local)) {
            type = get_pair_type(idx[jj - 1] + k, ptype);

            mm3 = S1[jj];
            en  = my_c[idx[jj - 1] + k];
            if (sn[k] != sn[jj - 1])
              en += P->DuplexInit;         /* ??? */

            if (sc)
              if (sc->energy_up)
                en += sc->energy_up[jj][1];

            if (fij == en + my_fc[k - 1] + E_ExtLoop(type, -1, mm3, P)) {
              bp_stack[++(*stack_count)].i  = k;
              bp_stack[(*stack_count)].j    = jj - 1;
              *u                            = k - 1;
              *i                            = k;
              *j                            = jj - 1;
              return 1;
            }
          }

          if ((k > ii) && (sn[jj - 1] == sn[jj]) &&
              evaluate(ii, jj, k - 2, k, VRNA_DECOMP_EXT_EXT_STEM1, &hc_dat_local)) {
            type = get_pair_type(idx[jj - 1] + k, ptype);

            mm3 = S1[jj];
            en  = my_c[idx[jj - 1] + k];
            if (sn[k] != sn[jj - 1])
              en += P->DuplexInit;         /* ??? */

            mm5 = (sn[k - 1] == sn[k]) ? S1[k - 1] : -1;
            if (sc)
              if (sc->energy_up)
                en += sc->energy_up[k - 1][1];

            if (fij == my_fc[k - 2] + en + E_ExtLoop(type, mm5, mm3, P)) {
              bp_stack[++(*stack_count)].i  = k;
              bp_stack[(*stack_count)].j    = jj - 1;
              *u                            = k - 2;
              *i                            = k;
              *j                            = jj - 1;
              return 1;
            }
          }
        }
        break;
    }
  }

  return 0;
}


PRIVATE int
BT_mb_loop_split(vrna_fold_compound_t *vc,
                 int                  *i,
                 int                  *j,
                 int                  *k,
                 int                  *l,
                 int                  *component1,
                 int                  *component2,
                 vrna_bp_stack_t      *bp_stack,
                 int                  *stack_count)
{
  char                      *ptype;
  short                     *S1;
  int                       ij, ii, jj, fij, fi, u, en, *my_c, *my_fML, *my_ggg,
                            turn, *idx, with_gquad, dangle_model, *rtype, kk, cnt,
                            with_ud, type, type_2;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_ud_t                 *domains_up;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  P           = vc->params;
  md          = &(P->model_details);
  hc          = vc->hc;
  sc          = vc->sc;
  idx         = vc->jindx;
  ptype       = vc->ptype;
  rtype       = &(md->rtype[0]);
  S1          = vc->sequence_encoding;
  domains_up  = vc->domains_up;

  my_c          = vc->matrices->c;
  my_fML        = vc->matrices->fML;
  my_ggg        = vc->matrices->ggg;
  turn          = md->min_loop_size;
  with_gquad    = md->gquad;
  with_ud       = (domains_up && domains_up->energy_cb) ? 1 : 0;
  dangle_model  = md->dangles;

  hc_dat_local.idx    = vc->jindx;
  hc_dat_local.mx     = hc->matrix;
  hc_dat_local.hc_up  = hc->up_ml;
  hc_dat_local.cp     = vc->cutpoint;

  if (hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default;
  }

  ii  = *i;
  jj  = *j;

  if (with_ud) {
    /* nibble off unpaired stretches at 3' site */
    do {
      fij = my_fML[idx[jj] + ii];
      fi  = INF;

      /* process regular unpaired nucleotides (unbound by ligand) first */
      if (evaluate(ii, jj, ii, jj - 1, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
        fi = my_fML[idx[jj - 1] + ii] +
             P->MLbase;

        if (sc) {
          if (sc->energy_up)
            fi += sc->energy_up[jj][1];

          if (sc->f)
            fi += sc->f(ii, jj, ii, jj - 1, VRNA_DECOMP_ML_ML, sc->data);
        }

        if (jj == ii)
          return 0; /* no more pairs */

        if (fij == fi) {
          jj--;
          continue;
        }
      }

      /* next try to nibble off ligand */
      for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
        u   = domains_up->uniq_motif_size[cnt];
        kk  = jj - u + 1;
        if ((kk >= ii) && evaluate(ii, jj, ii, jj - u, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
          en = domains_up->energy_cb(vc,
                                     kk,
                                     jj,
                                     VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                     domains_up->data);

          if (sc) {
            if (sc->energy_up)
              en += sc->energy_up[kk][u];

            if (sc->f)
              en += sc->f(ii, jj, ii, jj - u, VRNA_DECOMP_ML_ML, sc->data);
          }

          fi = my_fML[idx[kk - 1] + ii] +
               u * P->MLbase;

          fi += en;

          if (fij == fi) {
            /* skip remaining motifs after first hit */
            jj = kk - 1;
            break;
          }
        }
      }

      if (jj < ii)
        return 0; /* no more pairs */
    } while (fij == fi);

    /* nibble off unpaired stretches at 5' site */
    do {
      fij = my_fML[idx[jj] + ii];
      fi  = INF;

      /* again, process regular unpaired nucleotides (unbound by ligand) first */
      if (evaluate(ii, jj, ii + 1, jj, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
        fi = my_fML[idx[jj] + ii + 1] +
             P->MLbase;

        if (sc) {
          if (sc->energy_up)
            fi += sc->energy_up[ii][1];

          if (sc->f)
            fi += sc->f(ii, jj, ii + 1, jj, VRNA_DECOMP_ML_ML, sc->data);
        }

        if (ii + 1 == jj)
          return 0; /* no more pairs */

        if (fij == fi) {
          ii++;
          continue;
        }
      }

      /* next try to nibble off ligand again */
      for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
        u   = domains_up->uniq_motif_size[cnt];
        kk  = ii + u - 1;
        if ((kk <= jj) && evaluate(ii, jj, ii + u, jj, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
          en = domains_up->energy_cb(vc,
                                     ii,
                                     kk,
                                     VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                     domains_up->data);

          if (sc) {
            if (sc->energy_up)
              en += sc->energy_up[ii][u];

            if (sc->f)
              en += sc->f(ii, jj, ii + u, jj, VRNA_DECOMP_ML_ML, sc->data);
          }

          fi = my_fML[idx[jj] + kk + 1] +
               u * P->MLbase;
          fi += en;

          if (fij == fi) {
            /* skip remaining motifs after first hit */
            ii = kk + 1;
            break;
          }
        }
      }

      if (ii > jj)
        return 0; /* no more pairs */
    } while (fij == fi);
  } else {
    /* nibble off unpaired 3' bases */
    do {
      fij = my_fML[idx[jj] + ii];
      fi  = INF;

      if (evaluate(ii, jj, ii, jj - 1, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
        fi = my_fML[idx[jj - 1] + ii] +
             P->MLbase;

        if (sc) {
          if (sc->energy_up)
            fi += sc->energy_up[jj][1];

          if (sc->f)
            fi += sc->f(ii, jj, ii, jj - 1, VRNA_DECOMP_ML_ML, sc->data);
        }
      }

      if (--jj == 0)
        break;
    } while (fij == fi);
    jj++;

    /* nibble off unpaired 5' bases */
    do {
      fij = my_fML[idx[jj] + ii];
      fi  = INF;

      if (evaluate(ii, jj, ii + 1, jj, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
        fi = my_fML[idx[jj] + ii + 1] +
             P->MLbase;

        if (sc) {
          if (sc->energy_up)
            fi += sc->energy_up[ii][1];

          if (sc->f)
            fi += sc->f(ii, jj, ii + 1, jj, VRNA_DECOMP_ML_ML, sc->data);
        }
      }

      if (++ii == jj)
        break;
    } while (fij == fi);
    ii--;

    if (jj < ii + turn + 1) /* no more pairs */
      return 0;
  }

  ij = idx[jj] + ii;

  *component1 = *component2 = 1; /* split into two multi loop parts by default */

  /* 1. test for single component */

  if (with_gquad) {
    if (fij == my_ggg[ij] + E_MLstem(0, -1, -1, P)) {
      *i  = *j = -1;
      *k  = *l = -1;
      vrna_BT_gquad_mfe(vc, ii, jj, bp_stack, stack_count);
      return 1;
    }
  }

  type  = get_pair_type(ij, ptype);
  en    = my_c[ij];

  if (sc)
    if (sc->f)
      en += sc->f(ii, jj, ii, jj, VRNA_DECOMP_ML_STEM, sc->data);

  switch (dangle_model) {
    case 0:
      if (evaluate(ii, jj, ii, jj, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        if (fij == en + E_MLstem(type, -1, -1, P)) {
          *i          = *j = -1;
          *k          = ii;
          *l          = jj;
          *component2 = 2;          /* 2nd part is structure enclosed by base pair */
          return 1;
        }
      }

      break;

    case 2:
      if (evaluate(ii, jj, ii, jj, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        if (fij == en + E_MLstem(type, S1[ii - 1], S1[jj + 1], P)) {
          *i          = *j = -1;
          *k          = ii;
          *l          = jj;
          *component2 = 2;
          return 1;
        }
      }

      break;

    default:
      if (evaluate(ii, jj, ii, jj, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        if (fij == en + E_MLstem(type, -1, -1, P)) {
          *i          = *j = -1;
          *k          = ii;
          *l          = jj;
          *component2 = 2;
          return 1;
        }
      }

      if (evaluate(ii, jj, ii + 1, jj, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        int tmp_en = fij;
        if (sc) {
          if (sc->energy_up)
            tmp_en -= sc->energy_up[ii][1];

          if (sc->f)
            tmp_en -= sc->f(ii, jj, ii + 1, jj, VRNA_DECOMP_ML_STEM, sc->data);
        }

        type = get_pair_type(ij + 1, ptype);

        if (tmp_en == my_c[ij + 1] + E_MLstem(type, S1[ii], -1, P) + P->MLbase) {
          *i          = *j = -1;
          *k          = ii + 1;
          *l          = jj;
          *component2 = 2;
          return 1;
        }
      }

      if (evaluate(ii, jj, ii, jj - 1, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        int tmp_en = fij;
        if (sc) {
          if (sc->energy_up)
            tmp_en -= sc->energy_up[jj][1];

          if (sc->f)
            tmp_en -= sc->f(ii, jj, ii, jj - 1, VRNA_DECOMP_ML_STEM, sc->data);
        }

        type = get_pair_type(idx[jj - 1] + ii, ptype);

        if (tmp_en == my_c[idx[jj - 1] + ii] + E_MLstem(type, -1, S1[jj], P) + P->MLbase) {
          *i          = *j = -1;
          *k          = ii;
          *l          = jj - 1;
          *component2 = 2;
          return 1;
        }
      }

      if (evaluate(ii, jj, ii + 1, jj - 1, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        int tmp_en = fij;
        if (sc) {
          if (sc->energy_up)
            tmp_en -= sc->energy_up[ii][1] + sc->energy_up[jj][1];

          if (sc->f)
            tmp_en -= sc->f(ii, jj, ii + 1, jj - 1, VRNA_DECOMP_ML_STEM, sc->data);
        }

        type = get_pair_type(idx[jj - 1] + ii + 1, ptype);

        if (tmp_en ==
            my_c[idx[jj - 1] + ii + 1] + E_MLstem(type, S1[ii], S1[jj], P) + 2 * P->MLbase) {
          *i          = *j = -1;
          *k          = ii + 1;
          *l          = jj - 1;
          *component2 = 2;
          return 1;
        }
      }

      break;
  }

  /* 2. Test for possible split point */
  if (hc->f) {
    for (u = ii + 1 + turn; u <= jj - 2 - turn; u++) {
      if (hc->f(ii, jj, u, u + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
        en = my_fML[idx[u] + ii] + my_fML[idx[jj] + u + 1];
        if (sc)
          if (sc->f)
            en += sc->f(ii, jj, u, u + 1, VRNA_DECOMP_ML_ML_ML, sc->data);

        if (fij == en) {
          *i  = ii;
          *j  = u;
          *k  = u + 1;
          *l  = jj;
          return 1;
        }
      }
    }
  } else {
    for (u = ii + 1 + turn; u <= jj - 2 - turn; u++) {
      en = my_fML[idx[u] + ii] + my_fML[idx[jj] + u + 1];
      if (sc)
        if (sc->f)
          en += sc->f(ii, jj, u, u + 1, VRNA_DECOMP_ML_ML_ML, sc->data);

      if (fij == en) {
        *i  = ii;
        *j  = u;
        *k  = u + 1;
        *l  = jj;
        return 1;
      }
    }
  }

  /* 3. last chance! Maybe coax stack */
  if (dangle_model == 3) {
    int ik, k1j, tmp_en;
    for (k1j = idx[jj] + ii + turn + 2, u = ii + 1 + turn; u <= jj - 2 - turn; u++, k1j++) {
      ik = idx[u] + ii;
      if (evaluate(ii, u, u + 1, jj, VRNA_DECOMP_ML_COAXIAL_ENC, &hc_dat_local)) {
        type    = rtype[get_pair_type(ik, ptype)];
        type_2  = rtype[get_pair_type(k1j, ptype)];

        tmp_en = my_c[ik] +
                 my_c[k1j] +
                 P->stack[type][type_2] +
                 2 * P->MLintern[1];

        if (sc)
          if (sc->f)
            tmp_en += sc->f(ii, u, u + 1, jj, VRNA_DECOMP_ML_COAXIAL, sc->data);

        if (fij == tmp_en) {
          *i          = ii;
          *j          = u;
          *k          = u + 1;
          *l          = jj;
          *component1 = *component2 = 2;
          return 1;
        }
      }
    }
  }

  return 0;
}


PRIVATE int
BT_mb_loop_split_comparative(vrna_fold_compound_t *vc,
                             int                  *i,
                             int                  *j,
                             int                  *k,
                             int                  *l,
                             int                  *component1,
                             int                  *component2,
                             vrna_bp_stack_t      *bp_stack,
                             int                  *stack_count)
{
  unsigned int              **a2s;
  short                     **S, **S5, **S3;
  int                       ij, ii, jj, fij, fi, u, en, *my_c, *my_fML, *my_ggg,
                            turn, *idx, with_gquad, dangle_model, ss, n_seq, tt;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n_seq = vc->n_seq;
  S     = vc->S;
  S5    = vc->S5;
  S3    = vc->S3;
  a2s   = vc->a2s;
  P     = vc->params;
  md    = &(P->model_details);
  hc    = vc->hc;
  scs   = vc->scs;
  idx   = vc->jindx;

  my_c          = vc->matrices->c;
  my_fML        = vc->matrices->fML;
  my_ggg        = vc->matrices->ggg;
  turn          = md->min_loop_size;
  with_gquad    = md->gquad;
  dangle_model  = md->dangles;

  hc_dat_local.idx    = vc->jindx;
  hc_dat_local.mx     = hc->matrix;
  hc_dat_local.hc_up  = hc->up_ml;
  hc_dat_local.cp     = vc->cutpoint;

  if (hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default;
  }

  ii  = *i;
  jj  = *j;

  /* nibble off unpaired 3' bases */
  do {
    fij = my_fML[idx[jj] + ii];
    fi  = INF;

    if (evaluate(ii, jj, ii, jj - 1, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
      fi = my_fML[idx[jj - 1] + ii] +
           n_seq * P->MLbase;

      if (scs) {
        for (ss = 0; ss < n_seq; ss++) {
          if (scs[ss]) {
            if (scs[ss]->energy_up)
              fi += scs[ss]->energy_up[a2s[ss][jj]][1];

            if (scs[ss]->f)
              fi += scs[ss]->f(ii, jj, ii, jj - 1, VRNA_DECOMP_ML_ML, scs[ss]->data);
          }
        }
      }
    }

    if (--jj == 0)
      break;
  } while (fij == fi);
  jj++;

  /* nibble off unpaired 5' bases */
  do {
    fij = my_fML[idx[jj] + ii];
    fi  = INF;

    if (evaluate(ii, jj, ii + 1, jj, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
      fi = my_fML[idx[jj] + ii + 1] +
           n_seq * P->MLbase;

      if (scs) {
        for (ss = 0; ss < n_seq; ss++) {
          if (scs[ss]) {
            if (scs[ss]->energy_up)
              fi += scs[ss]->energy_up[a2s[ss][ii]][1];

            if (scs[ss]->f)
              fi += scs[ss]->f(ii, jj, ii + 1, jj, VRNA_DECOMP_ML_ML, scs[ss]->data);
          }
        }
      }
    }

    if (++ii == jj)
      break;
  } while (fij == fi);
  ii--;

  if (jj < ii + turn + 1) /* no more pairs */
    return 0;

  ij = idx[jj] + ii;

  *component1 = *component2 = 1; /* split into two multi loop parts by default */

  /* 1. test for single component */

  if (with_gquad) {
    if (fij == my_ggg[ij] + n_seq * E_MLstem(0, -1, -1, P)) {
      *i  = *j = -1;
      *k  = *l = -1;
      vrna_BT_gquad_mfe(vc, ii, jj, bp_stack, stack_count);
      return 1;
    }
  }

  en = my_c[ij];

  if (scs) {
    for (ss = 0; ss < n_seq; ss++) {
      if (scs[ss])
        if (scs[ss]->f)
          en += scs[ss]->f(ii, jj, ii, jj, VRNA_DECOMP_ML_STEM, scs[ss]->data);
    }
  }

  switch (dangle_model) {
    case 0:
      if (evaluate(ii, jj, ii, jj, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        for (ss = 0; ss < n_seq; ss++) {
          tt  = get_pair_type_md(S[ss][ii], S[ss][jj], md);
          en  += E_MLstem(tt, -1, -1, P);
        }
      }

      break;

    case 2:
      if (evaluate(ii, jj, ii, jj, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        for (ss = 0; ss < n_seq; ss++) {
          tt  = get_pair_type_md(S[ss][ii], S[ss][jj], md);
          en  += E_MLstem(tt, S5[ss][ii], S3[ss][jj], P);
        }
      }

      break;
  }

  if (fij == en) {
    *i          = *j = -1;
    *k          = ii;
    *l          = jj;
    *component2 = 2;          /* 2nd part is structure enclosed by base pair */
    return 1;
  }

  /* 2. Test for possible split point */
  if (hc->f) {
    for (u = ii + 1 + turn; u <= jj - 2 - turn; u++) {
      if (hc->f(ii, jj, u, u + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
        en = my_fML[idx[u] + ii] +
             my_fML[idx[jj] + u + 1];

        if (scs) {
          for (ss = 0; ss < n_seq; ss++)
            if (scs[ss] && scs[ss]->f)
              en += scs[ss]->f(ii, jj, u, u + 1, VRNA_DECOMP_ML_ML_ML, scs[ss]->data);
        }

        if (fij == en) {
          *i  = ii;
          *j  = u;
          *k  = u + 1;
          *l  = jj;
          return 1;
        }
      }
    }
  } else {
    for (u = ii + 1 + turn; u <= jj - 2 - turn; u++) {
      en = my_fML[idx[u] + ii] +
           my_fML[idx[jj] + u + 1];

      if (scs) {
        for (ss = 0; ss < n_seq; ss++)
          if (scs[ss] && scs[ss]->f)
            en += scs[ss]->f(ii, jj, u, u + 1, VRNA_DECOMP_ML_ML_ML, scs[ss]->data);
      }

      if (fij == en) {
        *i  = ii;
        *j  = u;
        *k  = u + 1;
        *l  = jj;
        return 1;
      }
    }
  }

  return 0;
}


PRIVATE int
BT_mb_loop_split_window(vrna_fold_compound_t  *vc,
                        int                   *i,
                        int                   *j,
                        int                   *k,
                        int                   *l,
                        int                   *component1,
                        int                   *component2,
                        vrna_bp_stack_t       *bp_stack,
                        int                   *stack_count)
{
  char                      **ptype;
  short                     *S1;
  int                       ii, jj, fij, fi, u, en, **c, **fML, **ggg, turn,
                            with_gquad, dangle_model, *rtype, type, type_2;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  P     = vc->params;
  md    = &(P->model_details);
  hc    = vc->hc;
  sc    = vc->sc;
  ptype = vc->ptype_local;
  rtype = &(md->rtype[0]);
  S1    = vc->sequence_encoding;

  c             = vc->matrices->c_local;
  fML           = vc->matrices->fML_local;
  ggg           = vc->matrices->ggg_local;
  turn          = md->min_loop_size;
  with_gquad    = md->gquad;
  dangle_model  = md->dangles;

  hc_dat_local.idx        = vc->jindx;
  hc_dat_local.mx         = hc->matrix;
  hc_dat_local.mx_window  = hc->matrix_local;
  hc_dat_local.hc_up      = hc->up_ml;
  hc_dat_local.cp         = vc->cutpoint;

  if (hc->f) {
    evaluate            = &hc_default_user_window;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default_window;
  }

  ii  = *i;
  jj  = *j;

  /* nibble off unpaired 3' bases */
  do {
    fij = fML[ii][jj - ii];
    fi  = INF;

    if (evaluate(ii, jj, ii, jj - 1, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
      fi = fML[ii][jj - 1 - ii] +
           P->MLbase;

      if (sc) {
        if (sc->energy_up)
          fi += sc->energy_up[jj][1];

        if (sc->f)
          fi += sc->f(ii, jj, ii, jj - 1, VRNA_DECOMP_ML_ML, sc->data);
      }
    }

    if (--jj == 0)
      break;
  } while (fij == fi);
  jj++;

  /* nibble off unpaired 5' bases */
  do {
    fij = fML[ii][jj - ii];
    fi  = INF;

    if (evaluate(ii, jj, ii + 1, jj, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
      fi = fML[ii + 1][jj - (ii + 1)] +
           P->MLbase;

      if (sc) {
        if (sc->energy_up)
          fi += sc->energy_up[ii][1];

        if (sc->f)
          fi += sc->f(ii, jj, ii + 1, jj, VRNA_DECOMP_ML_ML, sc->data);
      }
    }

    if (++ii == jj)
      break;
  } while (fij == fi);
  ii--;

  if (jj < ii + turn + 1) /* no more pairs */
    return 0;

  *component1 = *component2 = 1; /* split into two multi loop parts by default */

  /* 1. test for single component */

  if (with_gquad) {
    if (fij == ggg[ii][jj - ii] + E_MLstem(0, -1, -1, P)) {
      *i  = *j = -1;
      *k  = *l = -1;
      vrna_BT_gquad_mfe(vc, ii, jj, bp_stack, stack_count);
      return 1;
    }
  }

  type  = get_pair_type_window(ii, jj, ptype);
  en    = c[ii][jj - ii];

  if (sc)
    if (sc->f)
      en += sc->f(ii, jj, ii, jj, VRNA_DECOMP_ML_STEM, sc->data);

  switch (dangle_model) {
    case 0:
      if (evaluate(ii, jj, ii, jj, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        if (fij == en + E_MLstem(type, -1, -1, P)) {
          *i          = *j = -1;
          *k          = ii;
          *l          = jj;
          *component2 = 2;          /* 2nd part is structure enclosed by base pair */
          return 1;
        }
      }

      break;

    case 2:
      if (evaluate(ii, jj, ii, jj, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        if (fij == en + E_MLstem(type, S1[ii - 1], S1[jj + 1], P)) {
          *i          = *j = -1;
          *k          = ii;
          *l          = jj;
          *component2 = 2;
          return 1;
        }
      }

      break;

    default:
      if (evaluate(ii, jj, ii, jj, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        if (fij == en + E_MLstem(type, -1, -1, P)) {
          *i          = *j = -1;
          *k          = ii;
          *l          = jj;
          *component2 = 2;
          return 1;
        }
      }

      if (evaluate(ii, jj, ii + 1, jj, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        int tmp_en = fij;
        if (sc) {
          if (sc->energy_up)
            tmp_en -= sc->energy_up[ii][1];

          if (sc->f)
            tmp_en -= sc->f(ii, jj, ii + 1, jj, VRNA_DECOMP_ML_STEM, sc->data);
        }

        type = get_pair_type_window(ii + 1, jj, ptype);

        if (tmp_en == c[ii + 1][jj - (ii + 1)] + E_MLstem(type, S1[ii], -1, P) + P->MLbase) {
          *i          = *j = -1;
          *k          = ii + 1;
          *l          = jj;
          *component2 = 2;
          return 1;
        }
      }

      if (evaluate(ii, jj, ii, jj - 1, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        int tmp_en = fij;
        if (sc) {
          if (sc->energy_up)
            tmp_en -= sc->energy_up[jj][1];

          if (sc->f)
            tmp_en -= sc->f(ii, jj, ii, jj - 1, VRNA_DECOMP_ML_STEM, sc->data);
        }

        type = get_pair_type_window(ii, jj - 1, ptype);

        if (tmp_en == c[ii][jj - 1 - ii] + E_MLstem(type, -1, S1[jj], P) + P->MLbase) {
          *i          = *j = -1;
          *k          = ii;
          *l          = jj - 1;
          *component2 = 2;
          return 1;
        }
      }

      if (evaluate(ii, jj, ii + 1, jj - 1, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        int tmp_en = fij;
        if (sc) {
          if (sc->energy_up)
            tmp_en -= sc->energy_up[ii][1] + sc->energy_up[jj][1];

          if (sc->f)
            tmp_en -= sc->f(ii, jj, ii + 1, jj - 1, VRNA_DECOMP_ML_STEM, sc->data);
        }

        type = get_pair_type_window(ii + 1, jj - 1, ptype);

        if (tmp_en ==
            c[ii + 1][jj - 1 - (ii + 1)] + E_MLstem(type, S1[ii], S1[jj], P) + 2 * P->MLbase) {
          *i          = *j = -1;
          *k          = ii + 1;
          *l          = jj - 1;
          *component2 = 2;
          return 1;
        }
      }

      break;
  }

  /* 2. Test for possible split point */
  if (hc->f) {
    for (u = ii + 1 + turn; u <= jj - 2 - turn; u++) {
      if (hc->f(ii, jj, u, u + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
        en = fML[ii][u - ii] + fML[u + 1][jj - (u + 1)];
        if (sc)
          if (sc->f)
            en += sc->f(ii, jj, u, u + 1, VRNA_DECOMP_ML_ML_ML, sc->data);

        if (fij == en) {
          *i  = ii;
          *j  = u;
          *k  = u + 1;
          *l  = jj;
          return 1;
        }
      }
    }
  } else {
    for (u = ii + 1 + turn; u <= jj - 2 - turn; u++) {
      en = fML[ii][u - ii] + fML[u + 1][jj - (u + 1)];
      if (sc)
        if (sc->f)
          en += sc->f(ii, jj, u, u + 1, VRNA_DECOMP_ML_ML_ML, sc->data);

      if (fij == en) {
        *i  = ii;
        *j  = u;
        *k  = u + 1;
        *l  = jj;
        return 1;
      }
    }
  }

  /* 3. last chance! Maybe coax stack */
  if (dangle_model == 3) {
    int tmp_en;
    for (u = ii + 1 + turn; u <= jj - 2 - turn; u++) {
      if (evaluate(ii, u, u + 1, jj, VRNA_DECOMP_ML_COAXIAL_ENC, &hc_dat_local)) {
        type    = rtype[get_pair_type_window(ii, u, ptype)];
        type_2  = rtype[get_pair_type_window(u + 1, jj, ptype)];

        tmp_en = c[ii][u - ii] +
                 c[u + 1][jj - (u + 1)] +
                 P->stack[type][type_2] +
                 2 * P->MLintern[1];

        if (sc)
          if (sc->f)
            tmp_en += sc->f(ii, u, u + 1, jj, VRNA_DECOMP_ML_COAXIAL, sc->data);

        if (fij == tmp_en) {
          *i          = ii;
          *j          = u;
          *k          = u + 1;
          *l          = jj;
          *component1 = *component2 = 2;
          return 1;
        }
      }
    }
  }

  return 0;
}


PRIVATE int
BT_mb_loop_split_window_comparative(vrna_fold_compound_t  *vc,
                                    int                   *i,
                                    int                   *j,
                                    int                   *k,
                                    int                   *l,
                                    int                   *component1,
                                    int                   *component2,
                                    vrna_bp_stack_t       *bp_stack,
                                    int                   *stack_count)
{
  short                     **S, **S5, **S3;
  int                       ii, jj, fij, fi, u, en, **c, **fML, **ggg, n_seq,
                            turn, with_gquad, dangle_model, cc, ss, length, type;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  length  = vc->length;
  n_seq   = vc->n_seq;
  S       = vc->S;
  S5      = vc->S5;   /* S5[s][i] holds next base 5' of i in sequence s */
  S3      = vc->S3;   /* Sl[s][i] holds next base 3' of i in sequence s */
  P       = vc->params;
  md      = &(P->model_details);
  hc      = vc->hc;
  scs     = vc->scs;

  c             = vc->matrices->c_local;
  fML           = vc->matrices->fML_local;
  ggg           = vc->matrices->ggg_local;
  turn          = md->min_loop_size;
  with_gquad    = md->gquad;
  dangle_model  = md->dangles;

  hc_dat_local.mx         = hc->matrix;
  hc_dat_local.mx_window  = hc->matrix_local;
  hc_dat_local.hc_up      = hc->up_ml;
  hc_dat_local.cp         = vc->cutpoint;

  if (hc->f) {
    evaluate            = &hc_default_user_window;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default_window;
  }

  ii  = *i;
  jj  = *j;

  /* nibble off unpaired 3' bases */
  do {
    fij = fML[ii][jj - ii];
    fi  = INF;

    if (evaluate(ii, jj, ii, jj - 1, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
      fi = fML[ii][jj - 1 - ii] +
           n_seq * P->MLbase;

      if (scs) {
        for (ss = 0; ss < n_seq; ss++)
          if (scs[ss]) {
            if (scs[ss]->energy_up)
              fi += scs[ss]->energy_up[jj][1];

            if (scs[ss]->f)
              fi += scs[ss]->f(ii, jj, ii, jj - 1, VRNA_DECOMP_ML_ML, scs[ss]->data);
          }
      }
    }

    if (--jj == 0)
      break;
  } while (fij == fi);
  jj++;

  /* nibble off unpaired 5' bases */
  do {
    fij = fML[ii][jj - ii];
    fi  = INF;

    if (evaluate(ii, jj, ii + 1, jj, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
      fi = fML[ii + 1][jj - (ii + 1)] +
           n_seq * P->MLbase;

      if (scs) {
        for (ss = 0; ss < n_seq; ss++)
          if (scs[ss]) {
            if (scs[ss]->energy_up)
              fi += scs[ss]->energy_up[ii][1];

            if (scs[ss]->f)
              fi += scs[ss]->f(ii, jj, ii + 1, jj, VRNA_DECOMP_ML_ML, scs[ss]->data);
          }
      }
    }

    if (++ii == jj)
      break;
  } while (fij == fi);
  ii--;

  if (jj < ii + turn + 1) /* no more pairs */
    return 0;

  *component1 = *component2 = 1; /* split into two multi loop parts by default */

  /* 1. test for single component */

  if (with_gquad) {
    if (fij == ggg[ii][jj - ii] + n_seq * E_MLstem(0, -1, -1, P)) {
      *i  = *j = -1;
      *k  = *l = -1;
      vrna_BT_gquad_mfe(vc, ii, jj, bp_stack, stack_count);
      return 1;
    }
  }

  en = c[ii][jj - ii];

  if (scs) {
    for (ss = 0; ss < n_seq; ss++)
      if (scs[ss])
        if (scs[ss]->f)
          en += scs[ss]->f(ii, jj, ii, jj, VRNA_DECOMP_ML_STEM, scs[ss]->data);
  }

  switch (dangle_model) {
    case 0:
      if (evaluate(ii, jj, ii, jj, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        cc = 0;
        for (ss = 0; ss < n_seq; ss++) {
          type  = get_pair_type_md(S[ss][ii], S[ss][jj], md);
          cc    += E_MLstem(type, -1, -1, P);
        }

        if (fij == en + cc) {
          *i          = *j = -1;
          *k          = ii;
          *l          = jj;
          *component2 = 2;          /* 2nd part is structure enclosed by base pair */
          return 1;
        }
      }

      break;

    case 2:
      if (evaluate(ii, jj, ii, jj, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
        cc = 0;
        for (ss = 0; ss < n_seq; ss++) {
          type  = get_pair_type_md(S[ss][ii], S[ss][jj], md);
          cc    += E_MLstem(type, (ii > 1) ? S5[ss][ii] : -1, (jj < length) ? S3[ss][jj] : -1, P);
        }

        if (fij == en + cc) {
          *i          = *j = -1;
          *k          = ii;
          *l          = jj;
          *component2 = 2;
          return 1;
        }
      }

      break;
  }

  /* 2. Test for possible split point */
  if (hc->f) {
    for (u = ii + 1 + turn; u <= jj - 2 - turn; u++) {
      if (hc->f(ii, jj, u, u + 1, VRNA_DECOMP_ML_ML_ML, hc->data)) {
        en = fML[ii][u - ii] + fML[u + 1][jj - (u + 1)];
        if (scs) {
          for (ss = 0; ss < n_seq; ss++)
            if (scs[ss])
              if (scs[ss]->f)
                en += scs[ss]->f(ii, jj, u, u + 1, VRNA_DECOMP_ML_ML_ML, scs[ss]->data);
        }

        if (fij == en) {
          *i  = ii;
          *j  = u;
          *k  = u + 1;
          *l  = jj;
          return 1;
        }
      }
    }
  } else {
    for (u = ii + 1 + turn; u <= jj - 2 - turn; u++) {
      en = fML[ii][u - ii] + fML[u + 1][jj - (u + 1)];
      if (scs) {
        for (ss = 0; ss < n_seq; ss++)
          if (scs[ss])
            if (scs[ss]->f)
              en += scs[ss]->f(ii, jj, u, u + 1, VRNA_DECOMP_ML_ML_ML, scs[ss]->data);
      }

      if (fij == en) {
        *i  = ii;
        *j  = u;
        *k  = u + 1;
        *l  = jj;
        return 1;
      }
    }
  }

  return 0;
}


PRIVATE int
BT_mb_loop(vrna_fold_compound_t *vc,
           int                  *i,
           int                  *j,
           int                  *k,
           int                  en,
           int                  *component1,
           int                  *component2)
{
  char                      *ptype;
  short                     s5, s3, *S1;
  unsigned int              *sn;
  int                       ij, p, q, r, e, tmp_en, cp, *idx, turn, dangle_model,
                            *my_c, *my_fML, *my_fc, *rtype, type, type_2, tt;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;
  vrna_callback_hc_evaluate *evaluate_ext;
  struct default_data       hc_dat_local_ext;

  cp            = vc->cutpoint;
  idx           = vc->jindx;
  ij            = idx[*j] + *i;
  S1            = vc->sequence_encoding;
  P             = vc->params;
  md            = &(P->model_details);
  sn            = vc->strand_number;
  hc            = vc->hc;
  sc            = vc->sc;
  my_c          = vc->matrices->c;
  my_fML        = vc->matrices->fML;
  my_fc         = vc->matrices->fc;
  turn          = md->min_loop_size;
  ptype         = vc->ptype;
  rtype         = &(md->rtype[0]);
  type          = get_pair_type(ij, ptype);
  tt            = type;
  type          = rtype[type];
  dangle_model  = md->dangles;

  hc_dat_local.idx    = vc->jindx;
  hc_dat_local.mx     = hc->matrix;
  hc_dat_local.hc_up  = hc->up_ml;
  hc_dat_local.cp     = vc->cutpoint;

  if (hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default;
  }

  hc_dat_local_ext.idx    = vc->jindx;
  hc_dat_local_ext.mx     = hc->matrix;
  hc_dat_local_ext.hc_up  = hc->up_ext;
  hc_dat_local_ext.cp     = vc->cutpoint;

  if (hc->f) {
    evaluate_ext            = &hc_default_user_ext;
    hc_dat_local_ext.hc_f   = hc->f;
    hc_dat_local_ext.hc_dat = hc->data;
  } else {
    evaluate_ext = &hc_default_ext;
  }

  p = *i + 1;
  q = *j - 1;

  r = q - turn - 1;

  if (evaluate(*i, *j, p, q, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
    /* is it a fake multi-loop? */
    /* NOTE: do we really want to evaluate it hard-constraint-wise as a multibranch loop? */
    if (sn[*i] != sn[*j]) {
      int ii, jj;
      ii  = jj = 0;
      e   = my_fc[p] + my_fc[q];
      if (sc)
        if (sc->energy_bp)
          e += sc->energy_bp[ij];

      s5  = (sn[q] == sn[*j]) ? S1[q] : -1;
      s3  = (sn[*i] == sn[p]) ? S1[p] : -1;

      switch (dangle_model) {
        case 0:
          if (en == e + E_ExtLoop(type, -1, -1, P))
            ii = p, jj = q;

          break;

        case 2:
          if (en == e + E_ExtLoop(type, s5, s3, P))
            ii = p, jj = q;

          break;

        default:
          if (en == e + E_ExtLoop(type, -1, -1, P)) {
            ii = p, jj = q;
            break;
          }

          if (evaluate_ext(p, q, p + 1, q, VRNA_DECOMP_EXT_EXT, &hc_dat_local_ext)) {
            e = my_fc[p + 1] + my_fc[q];
            if (sc) {
              if (sc->energy_up)
                e += sc->energy_up[p][1];

              if (sc->energy_bp)
                e += sc->energy_bp[ij];
            }

            if (en == e + E_ExtLoop(type, -1, s3, P)) {
              ii  = p + 1;
              jj  = q;
              break;
            }
          }

          if (evaluate_ext(p, q, p, q - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local_ext)) {
            e = my_fc[p] + my_fc[q - 1];
            if (sc) {
              if (sc->energy_up)
                e += sc->energy_up[q][1];

              if (sc->energy_bp)
                e += sc->energy_bp[ij];
            }

            if (en == e + E_ExtLoop(type, s5, -1, P)) {
              ii  = p;
              jj  = q - 1;
              break;
            }
          }

          if (evaluate_ext(p, q, p + 1, q - 1, VRNA_DECOMP_EXT_EXT, &hc_dat_local_ext)) {
            e = my_fc[p + 1] + my_fc[q - 1];
            if (sc) {
              if (sc->energy_up)
                e += sc->energy_up[p][1] + sc->energy_up[q][1];

              if (sc->energy_bp)
                e += sc->energy_bp[ij];
            }

            if (en == e + E_ExtLoop(type, s5, s3, P)) {
              ii  = p + 1;
              jj  = q - 1;
              break;
            }
          }

          break;
      }

      if (ii) {
        /* found a decomposition */
        *component1 = 3;
        *i          = ii;
        *k          = cp - 1;
        *j          = jj;
        *component2 = 4;
        return 1;
      }
    }

    /* true multi loop? */
    *component1 = *component2 = 1;  /* both components are MB loop parts by default */

    s5  = (sn[q] == sn[*j]) ? S1[q] : -1;
    s3  = (sn[*i] == sn[p]) ? S1[p] : -1;

    switch (dangle_model) {
      case 0:
        e = en - E_MLstem(type, -1, -1, P) - P->MLclosing;
        if (sc) {
          if (sc->energy_bp)
            e -= sc->energy_bp[ij];

          if (sc->f)
            e -= sc->f(*i, *j, p, q, VRNA_DECOMP_PAIR_ML, sc->data);
        }

        for (r = *i + 2 + turn; r < *j - 2 - turn; ++r) {
          if (evaluate(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
            tmp_en = my_fML[idx[r] + p] + my_fML[idx[q] + r + 1];
            if (sc)
              if (sc->f)
                tmp_en += sc->f(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, sc->data);

            if (e == tmp_en)
              break;
          }
        }
        break;

      case 2:
        e = en - E_MLstem(type, s5, s3, P) - P->MLclosing;
        if (sc) {
          if (sc->energy_bp)
            e -= sc->energy_bp[ij];

          if (sc->f)
            e -= sc->f(*i, *j, p, q, VRNA_DECOMP_PAIR_ML, sc->data);
        }

        for (r = p + turn + 1; r < q - turn - 1; ++r) {
          if (evaluate(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
            tmp_en = my_fML[idx[r] + p] + my_fML[idx[q] + r + 1];
            if (sc)
              if (sc->f)
                tmp_en += sc->f(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, sc->data);

            if (e == tmp_en)
              break;
          }
        }
        break;

      default:
        e = en - P->MLclosing;
        if (sc)
          if (sc->energy_bp)
            e -= sc->energy_bp[ij];

        for (r = p + turn + 1; r < q - turn - 1; ++r) {
          if (evaluate(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
            tmp_en = my_fML[idx[r] + p] + my_fML[idx[q] + r + 1] + E_MLstem(type, -1, -1, P);
            if (sc) {
              if (sc->f) {
                tmp_en  += sc->f(*i, *j, p, q, VRNA_DECOMP_PAIR_ML, sc->data);
                tmp_en  += sc->f(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
              }
            }

            if (e == tmp_en)
              break;
          }

          if (evaluate(*i, *j, p + 1, q, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
            if (evaluate(p + 1, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
              tmp_en = e;
              if (sc) {
                if (sc->energy_up)
                  tmp_en -= sc->energy_up[p][1];

                if (sc->f) {
                  tmp_en  -= sc->f(*i, *j, p + 1, q, VRNA_DECOMP_PAIR_ML, sc->data);
                  tmp_en  -= sc->f(p + 1, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
                }
              }

              if (tmp_en ==
                  my_fML[idx[r] + p + 1] + my_fML[idx[q] + r + 1] +
                  E_MLstem(type, -1, s3, P) + P->MLbase) {
                p += 1;
                break;
              }
            }
          }

          if (evaluate(*i, *j, p, q - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
            if (evaluate(p, q - 1, r, r + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
              tmp_en = e;
              if (sc) {
                if (sc->energy_up)
                  tmp_en -= sc->energy_up[q][1];

                if (sc->f) {
                  tmp_en  -= sc->f(*i, *j, p, q - 1, VRNA_DECOMP_PAIR_ML, sc->data);
                  tmp_en  -= sc->f(p, q - 1, r, r + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
                }
              }

              if (tmp_en ==
                  my_fML[idx[r] + p] + my_fML[idx[q - 1] + r + 1] +
                  E_MLstem(type, s5, -1, P) + P->MLbase) {
                q -= 1;
                break;
              }
            }
          }

          if (evaluate(*i, *j, p + 1, q - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
            if (evaluate(p + 1, q - 1, r, r + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
              tmp_en = e;
              if (sc) {
                if (sc->energy_up)
                  tmp_en -= sc->energy_up[p][1] + sc->energy_up[q][1];

                if (sc->f) {
                  tmp_en  -= sc->f(*i, *j, p + 1, q - 1, VRNA_DECOMP_PAIR_ML, sc->data);
                  tmp_en  -= sc->f(p + 1, q - 1, r, r + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
                }
              }

              if (tmp_en ==
                  my_fML[idx[r] + p + 1] + my_fML[idx[q - 1] + r + 1] +
                  E_MLstem(type, s5, s3, P) + 2 * P->MLbase) {
                p += 1;
                q -= 1;
                break;
              }
            }
          }

          /* coaxial stacking of (i.j) with (i+1.r) or (r.j-1) */
          /* use MLintern[1] since coax stacked pairs don't get TerminalAU */
          if (dangle_model == 3) {
            tmp_en = e;
            if (evaluate(*i, *j, p, r, VRNA_DECOMP_ML_COAXIAL, &hc_dat_local)) {
              type_2 = rtype[get_pair_type(idx[r] + p, ptype)];

              tmp_en = my_c[idx[r] + p] +
                       P->stack[tt][type_2] +
                       my_fML[idx[q] + r + 1];

              if (sc) {
                if (sc->f) {
                  tmp_en  += sc->f(*i, *j, p, q, VRNA_DECOMP_PAIR_ML, sc->data);
                  tmp_en  += sc->f(*i, *j, p, r, VRNA_DECOMP_ML_COAXIAL, sc->data);
                }
              }

              if (e == tmp_en + 2 * P->MLintern[1]) {
                *component1 = 2;
                break;
              }
            }

            if (evaluate(*i, *j, r + 1, q, VRNA_DECOMP_ML_COAXIAL, &hc_dat_local)) {
              type_2 = rtype[get_pair_type(idx[q] + r + 1, ptype)];

              tmp_en = my_c[idx[q] + r + 1] +
                       P->stack[tt][type_2] +
                       my_fML[idx[r] + p];

              if (sc) {
                if (sc->f) {
                  tmp_en  += sc->f(*i, *j, p, q, VRNA_DECOMP_PAIR_ML, sc->data);
                  tmp_en  += sc->f(*i, *j, r + 1, q, VRNA_DECOMP_ML_COAXIAL, sc->data);
                }
              }

              if (e == tmp_en + 2 * P->MLintern[1]) {
                *component2 = 2;
                break;
              }
            }
          }
        }
        break;
    }
  }

  if (r <= *j - turn - 3) {
    *i  = p;
    *k  = r;
    *j  = q;
    return 1;
  } else {
#if 0
    /* Y shaped ML loops fon't work yet */
    if (dangle_model == 3) {
      d5  = P->dangle5[tt][S1[j - 1]];
      d3  = P->dangle3[tt][S1[i + 1]];
      /* (i,j) must close a Y shaped ML loop with coax stacking */
      if (cij == fML[indx[j - 2] + i + 2] + mm + d3 + d5 + P->MLbase + P->MLbase) {
        i1  = i + 2;
        j1  = j - 2;
      } else if (cij == fML[indx[j - 2] + i + 1] + mm + d5 + P->MLbase) {
        j1 = j - 2;
      } else if (cij == fML[indx[j - 1] + i + 2] + mm + d3 + P->MLbase) {
        i1 = i + 2;
      } else /* last chance */
      if (cij != fML[indx[j - 1] + i + 1] + mm + P->MLbase) {
        fprintf(stderr, "backtracking failed in repeat");
      }

      /* if we arrive here we can express cij via fML[i1,j1]+dangles */
      bt_stack[++s].i = i1;
      bt_stack[s].j   = j1;
    }

#endif
  }

  return 0;
}


PRIVATE int
BT_mb_loop_comparative(vrna_fold_compound_t *vc,
                       int                  *i,
                       int                  *j,
                       int                  *k,
                       int                  en,
                       int                  *component1,
                       int                  *component2)
{
  short                     **S, **S5, **S3;
  int                       ij, p, q, r, e, tmp_en, *idx, turn, dangle_model,
                            *my_fML, ss, n_seq, tt;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n_seq         = vc->n_seq;
  S             = vc->S;
  S5            = vc->S5;
  S3            = vc->S3;
  idx           = vc->jindx;
  ij            = idx[*j] + *i;
  P             = vc->params;
  md            = &(P->model_details);
  hc            = vc->hc;
  scs           = vc->scs;
  my_fML        = vc->matrices->fML;
  turn          = md->min_loop_size;
  dangle_model  = md->dangles;

  hc_dat_local.idx    = vc->jindx;
  hc_dat_local.mx     = hc->matrix;
  hc_dat_local.hc_up  = hc->up_ml;
  hc_dat_local.cp     = vc->cutpoint;

  if (hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default;
  }

  p = *i + 1;
  q = *j - 1;

  r = q - turn - 1;

  if (evaluate(*i, *j, p, q, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
    *component1 = *component2 = 1;  /* both components are MB loop parts by default */

    e = en - n_seq * P->MLclosing;

    switch (dangle_model) {
      case 0:
        for (ss = 0; ss < n_seq; ss++) {
          tt  = get_pair_type_md(S[ss][*j], S[ss][*i], md);
          e   -= E_MLstem(tt, -1, -1, P);
        }
        break;

      case 2:
        for (ss = 0; ss < n_seq; ss++) {
          tt  = get_pair_type_md(S[ss][*j], S[ss][*i], md);
          e   -= E_MLstem(tt, S5[ss][*j], S3[ss][*i], P);
        }
        break;
    }

    if (scs) {
      for (ss = 0; ss < n_seq; ss++) {
        if (scs[ss]) {
          if (scs[ss]->energy_bp)
            e -= scs[ss]->energy_bp[ij];

          if (scs[ss]->f)
            e -= scs[ss]->f(*i, *j, p, q, VRNA_DECOMP_PAIR_ML, scs[ss]->data);
        }
      }
    }

    for (r = p + turn + 1; r < q - turn - 1; ++r) {
      if (evaluate(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
        tmp_en = my_fML[idx[r] + p] + my_fML[idx[q] + r + 1];
        if (scs) {
          for (ss = 0; ss < n_seq; ss++) {
            if (scs[ss])
              if (scs[ss]->f)
                tmp_en += scs[ss]->f(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, scs[ss]->data);
          }
        }

        if (e == tmp_en)
          break;
      }
    }
  }

  if (r <= *j - turn - 3) {
    *i  = p;
    *k  = r;
    *j  = q;
    return 1;
  }

  return 0;
}


PRIVATE int
BT_mb_loop_window(vrna_fold_compound_t  *vc,
                  int                   *i,
                  int                   *j,
                  int                   *k,
                  int                   en,
                  int                   *component1,
                  int                   *component2)
{
  char                      **ptype;
  short                     s5, s3, *S1;
  unsigned int              *sn;
  int                       p, q, r, e, tmp_en, turn, dangle_model,
                            **c, **fML, *rtype, type, type_2, tt;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  S1            = vc->sequence_encoding;
  P             = vc->params;
  md            = &(P->model_details);
  sn            = vc->strand_number;
  hc            = vc->hc;
  sc            = vc->sc;
  c             = vc->matrices->c_local;
  fML           = vc->matrices->fML_local;
  turn          = md->min_loop_size;
  ptype         = vc->ptype_local;
  rtype         = &(md->rtype[0]);
  type          = get_pair_type_window(*i, *j, ptype);
  tt            = type;
  type          = rtype[type];
  dangle_model  = md->dangles;

  hc_dat_local.idx        = vc->jindx;
  hc_dat_local.mx         = hc->matrix;
  hc_dat_local.mx_window  = hc->matrix_local;
  hc_dat_local.hc_up      = hc->up_ml;
  hc_dat_local.cp         = vc->cutpoint;

  if (hc->f) {
    evaluate            = &hc_default_user_window;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default_window;
  }

  p = *i + 1;
  q = *j - 1;

  r = q - turn - 1;

  if (evaluate(*i, *j, p, q, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
    *component1 = *component2 = 1;  /* both components are MB loop parts by default */

    s5  = (sn[q] == sn[*j]) ? S1[q] : -1;
    s3  = (sn[*i] == sn[p]) ? S1[p] : -1;

    switch (dangle_model) {
      case 0:
        e = en - E_MLstem(type, -1, -1, P) - P->MLclosing;
        if (sc) {
          if (sc->energy_bp_local)
            e -= sc->energy_bp_local[*i][*j - *i];

          if (sc->f)
            e -= sc->f(*i, *j, p, q, VRNA_DECOMP_PAIR_ML, sc->data);
        }

        for (r = *i + 2 + turn; r < *j - 2 - turn; ++r) {
          if (evaluate(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
            tmp_en = fML[p][r - p] +
                     fML[r + 1][q - (r + 1)];

            if (sc)
              if (sc->f)
                tmp_en += sc->f(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, sc->data);

            if (e == tmp_en)
              break;
          }
        }
        break;

      case 2:
        e = en - E_MLstem(type, s5, s3, P) - P->MLclosing;
        if (sc) {
          if (sc->energy_bp_local)
            e -= sc->energy_bp_local[*i][*j - *i];

          if (sc->f)
            e -= sc->f(*i, *j, p, q, VRNA_DECOMP_PAIR_ML, sc->data);
        }

        for (r = p + turn + 1; r < q - turn - 1; ++r) {
          if (evaluate(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
            tmp_en = fML[p][r - p] +
                     fML[r + 1][q - (r + 1)];

            if (sc)
              if (sc->f)
                tmp_en += sc->f(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, sc->data);

            if (e == tmp_en)
              break;
          }
        }
        break;

      default:
        e = en - P->MLclosing;

        for (r = p + turn + 1; r < q - turn - 1; ++r) {
          if (evaluate(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
            tmp_en = fML[p][r - p] +
                     fML[r + 1][q - (r + 1)] +
                     E_MLstem(type, -1, -1, P);

            if (sc) {
              if (sc->energy_bp)
                tmp_en += sc->energy_bp_local[*i][*j - *i];

              if (sc->f) {
                tmp_en  += sc->f(*i, *j, p, q, VRNA_DECOMP_PAIR_ML, sc->data);
                tmp_en  += sc->f(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
              }
            }

            if (e == tmp_en)
              break;
          }

          if (evaluate(*i, *j, p + 1, q, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
            if (evaluate(p + 1, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
              tmp_en = fML[p + 1][r - (p + 1)] +
                       fML[r + 1][q - (r + 1)] +
                       E_MLstem(type, -1, s3, P) +
                       P->MLbase;

              if (sc) {
                if (sc->energy_up)
                  tmp_en += sc->energy_up[p][1];

                if (sc->energy_bp)
                  tmp_en += sc->energy_bp_local[*i][*j - *i];

                if (sc->f) {
                  tmp_en  += sc->f(*i, *j, p + 1, q, VRNA_DECOMP_PAIR_ML, sc->data);
                  tmp_en  += sc->f(p + 1, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
                }
              }

              if (e == tmp_en) {
                p += 1;
                break;
              }
            }
          }

          if (evaluate(*i, *j, p, q - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
            if (evaluate(p, q - 1, r, r + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
              tmp_en = fML[p][r - p] +
                       fML[r + 1][q - 1 - (r + 1)] +
                       E_MLstem(type, s5, -1, P) +
                       P->MLbase;

              if (sc) {
                if (sc->energy_up)
                  tmp_en += sc->energy_up[q][1];

                if (sc->energy_bp)
                  tmp_en += sc->energy_bp_local[*i][*j - *i];

                if (sc->f) {
                  tmp_en  += sc->f(*i, *j, p, q - 1, VRNA_DECOMP_PAIR_ML, sc->data);
                  tmp_en  += sc->f(p, q - 1, r, r + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
                }
              }

              if (e == tmp_en) {
                q -= 1;
                break;
              }
            }
          }

          if (evaluate(*i, *j, p + 1, q - 1, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
            if (evaluate(p + 1, q - 1, r, r + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
              tmp_en = fML[p + 1][r - (p + 1)] +
                       fML[r + 1][q - 1 - (r + 1)] +
                       E_MLstem(type, s5, s3, P) +
                       2 * P->MLbase;

              if (sc) {
                if (sc->energy_up)
                  tmp_en += sc->energy_up[p][1] +
                            sc->energy_up[q][1];

                if (sc->energy_bp)
                  tmp_en += sc->energy_bp_local[*i][*j - *i];

                if (sc->f) {
                  tmp_en  += sc->f(*i, *j, p + 1, q - 1, VRNA_DECOMP_PAIR_ML, sc->data);
                  tmp_en  += sc->f(p + 1, q - 1, r, r + 1, VRNA_DECOMP_ML_ML_ML, sc->data);
                }
              }

              if (e == tmp_en) {
                p += 1;
                q -= 1;
                break;
              }
            }
          }

          /* coaxial stacking of (i.j) with (i+1.r) or (r.j-1) */
          /* use MLintern[1] since coax stacked pairs don't get TerminalAU */
          if (dangle_model == 3) {
            tmp_en = e;
            if (evaluate(*i, *j, p, r, VRNA_DECOMP_ML_COAXIAL, &hc_dat_local)) {
              type_2 = rtype[get_pair_type_window(p, r, ptype)];

              tmp_en = c[p][r - p] +
                       P->stack[tt][type_2] +
                       fML[r + 1][q - (r + 1)];

              if (sc) {
                if (sc->energy_bp)
                  tmp_en += sc->energy_bp_local[*i][*j - *i];

                if (sc->f) {
                  tmp_en  += sc->f(*i, *j, p, q, VRNA_DECOMP_PAIR_ML, sc->data);
                  tmp_en  += sc->f(*i, *j, p, r, VRNA_DECOMP_ML_COAXIAL, sc->data);
                }
              }

              if (e == tmp_en + 2 * P->MLintern[1]) {
                *component1 = 2;
                break;
              }
            }

            if (evaluate(*i, *j, r + 1, q, VRNA_DECOMP_ML_COAXIAL, &hc_dat_local)) {
              type_2 = rtype[get_pair_type_window(r + 1, q, ptype)];

              tmp_en = c[r + 1][q - (r + 1)] +
                       P->stack[tt][type_2] +
                       fML[p][r - p];

              if (sc) {
                if (sc->energy_bp)
                  tmp_en += sc->energy_bp_local[*i][*j - *i];

                if (sc->f) {
                  tmp_en  += sc->f(*i, *j, p, q, VRNA_DECOMP_PAIR_ML, sc->data);
                  tmp_en  += sc->f(*i, *j, r + 1, q, VRNA_DECOMP_ML_COAXIAL, sc->data);
                }
              }

              if (e == tmp_en + 2 * P->MLintern[1]) {
                *component2 = 2;
                break;
              }
            }
          }
        }
        break;
    }
  }

  if (r <= *j - turn - 3) {
    *i  = p;
    *k  = r;
    *j  = q;
    return 1;
  } else {
#if 0
    /* Y shaped ML loops fon't work yet */
    if (dangle_model == 3) {
      d5  = P->dangle5[tt][S1[j - 1]];
      d3  = P->dangle3[tt][S1[i + 1]];
      /* (i,j) must close a Y shaped ML loop with coax stacking */
      if (cij == fML[indx[j - 2] + i + 2] + mm + d3 + d5 + P->MLbase + P->MLbase) {
        i1  = i + 2;
        j1  = j - 2;
      } else if (cij == fML[indx[j - 2] + i + 1] + mm + d5 + P->MLbase) {
        j1 = j - 2;
      } else if (cij == fML[indx[j - 1] + i + 2] + mm + d3 + P->MLbase) {
        i1 = i + 2;
      } else /* last chance */
      if (cij != fML[indx[j - 1] + i + 1] + mm + P->MLbase) {
        fprintf(stderr, "backtracking failed in repeat");
      }

      /* if we arrive here we can express cij via fML[i1,j1]+dangles */
      bt_stack[++s].i = i1;
      bt_stack[s].j   = j1;
    }

#endif
  }

  return 0;
}


PRIVATE int
BT_mb_loop_window_comparative(vrna_fold_compound_t  *vc,
                              int                   *i,
                              int                   *j,
                              int                   *k,
                              int                   en,
                              int                   *component1,
                              int                   *component2)
{
  short                     **S, **S5, **S3;
  int                       tt, p, q, r, e, tmp_en, turn, dangle_model,
                            **fML, mm, n_seq, ss;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n_seq         = vc->n_seq;
  S             = vc->S;
  S5            = vc->S5;   /* S5[s][i] holds next base 5' of i in sequence s */
  S3            = vc->S3;   /* Sl[s][i] holds next base 3' of i in sequence s */
  P             = vc->params;
  md            = &(P->model_details);
  hc            = vc->hc;
  scs           = vc->scs;
  fML           = vc->matrices->fML_local;
  turn          = md->min_loop_size;
  dangle_model  = md->dangles;

  hc_dat_local.idx        = vc->jindx;
  hc_dat_local.mx         = hc->matrix;
  hc_dat_local.mx_window  = hc->matrix_local;
  hc_dat_local.hc_up      = hc->up_ml;
  hc_dat_local.cp         = vc->cutpoint;

  if (hc->f) {
    evaluate            = &hc_default_user_window;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default_window;
  }

  p = *i + 1;
  q = *j - 1;

  r = q - turn - 1;

  if (evaluate(*i, *j, p, q, VRNA_DECOMP_PAIR_ML, &hc_dat_local)) {
    mm = n_seq * P->MLclosing;

    *component1 = *component2 = 1;  /* both components are MB loop parts by default */

    switch (dangle_model) {
      case 0:
        for (ss = 0; ss < n_seq; ss++) {
          tt  = get_pair_type_md(S[ss][*j], S[ss][*i], md);
          mm  += E_MLstem(tt, -1, -1, P);
        }

        e = en - mm;
        if (scs) {
          for (ss = 0; ss < n_seq; ss++)
            if (scs[ss])
              if (scs[ss]->f)
                e -= scs[ss]->f(*i, *j, p, q, VRNA_DECOMP_PAIR_ML, scs[ss]->data);
        }

        for (r = *i + 2 + turn; r < *j - 2 - turn; ++r) {
          if (evaluate(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
            tmp_en = fML[p][r - p] +
                     fML[r + 1][q - (r + 1)];

            if (scs) {
              for (ss = 0; ss < n_seq; ss++)
                if (scs[ss])
                  if (scs[ss]->f)
                    tmp_en += scs[ss]->f(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, scs[ss]->data);
            }

            if (e == tmp_en)
              break;
          }
        }
        break;

      case 2:
        for (ss = 0; ss < n_seq; ss++) {
          tt  = get_pair_type_md(S[ss][*j], S[ss][*i], md);
          mm  += E_MLstem(tt, S5[ss][*j], S3[ss][*i], P);
        }

        e = en - mm;
        if (scs) {
          for (ss = 0; ss < n_seq; ss++)
            if (scs[ss])
              if (scs[ss]->f)
                e -= scs[ss]->f(*i, *j, p, q, VRNA_DECOMP_PAIR_ML, scs[ss]->data);
        }

        for (r = p + turn + 1; r < q - turn - 1; ++r) {
          if (evaluate(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, &hc_dat_local)) {
            tmp_en = fML[p][r - p] +
                     fML[r + 1][q - (r + 1)];

            if (scs) {
              for (ss = 0; ss < n_seq; ss++)
                if (scs[ss])
                  if (scs[ss]->f)
                    tmp_en += scs[ss]->f(p, q, r, r + 1, VRNA_DECOMP_ML_ML_ML, scs[ss]->data);
            }

            if (e == tmp_en)
              break;
          }
        }
        break;
    }
  }

  if (r <= *j - turn - 3) {
    *i  = p;
    *k  = r;
    *j  = q;
    return 1;
  }

  return 0;
}


PUBLIC vrna_mx_pf_aux_ml_t *
vrna_exp_E_ml_fast_init(vrna_fold_compound_t *vc)
{
  vrna_mx_pf_aux_ml_t *aux_mx = NULL;

  if (vc) {
    int         i, j, d, n, u, turn, ij, *iidx;
    FLT_OR_DBL  *qm;

    n     = (int)vc->length;
    iidx  = vc->iindx;
    turn  = vc->exp_params->model_details.min_loop_size;
    qm    = vc->exp_matrices->qm;

    /* allocate memory for helper arrays */
    aux_mx            = (vrna_mx_pf_aux_ml_t *)vrna_alloc(sizeof(vrna_mx_pf_aux_ml_t));
    aux_mx->qqm       = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    aux_mx->qqm1      = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    aux_mx->qqmu_size = 0;
    aux_mx->qqmu      = NULL;

    if (vc->type == VRNA_FC_TYPE_SINGLE) {
      vrna_ud_t *domains_up = vc->domains_up;
      int       with_ud     = (domains_up && domains_up->exp_energy_cb);
      int       ud_max_size = 0;

      /* pre-processing ligand binding production rule(s) and auxiliary memory */
      if (with_ud) {
        for (u = 0; u < domains_up->uniq_motif_count; u++)
          if (ud_max_size < domains_up->uniq_motif_size[u])
            ud_max_size = domains_up->uniq_motif_size[u];

        aux_mx->qqmu_size = ud_max_size;
        aux_mx->qqmu      = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * (ud_max_size + 1));
        for (u = 0; u <= ud_max_size; u++)
          aux_mx->qqmu[u] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
      }
    }

    if (vc->hc->type == VRNA_HC_WINDOW) {
    } else {
      for (d = 0; d <= turn; d++)
        for (i = 1; i <= n - d; i++) {
          j   = i + d;
          ij  = iidx[i] - j;

          if (j > n)
            continue;

          qm[ij] = 0.;
        }
    }
  }

  return aux_mx;
}


PUBLIC void
vrna_exp_E_ml_fast_rotate(vrna_fold_compound_t  *vc,
                          vrna_mx_pf_aux_ml_t   *aux_mx)
{
  if (vc && aux_mx) {
    int         u;
    FLT_OR_DBL  *tmp;

    tmp           = aux_mx->qqm1;
    aux_mx->qqm1  = aux_mx->qqm;
    aux_mx->qqm   = tmp;

    /* rotate auxiliary arrays for unstructured domains */
    if (aux_mx->qqmu) {
      tmp = aux_mx->qqmu[aux_mx->qqmu_size];
      for (u = aux_mx->qqmu_size; u > 0; u--)
        aux_mx->qqmu[u] = aux_mx->qqmu[u - 1];
      aux_mx->qqmu[0] = tmp;
    }
  }
}


PUBLIC void
vrna_exp_E_ml_fast_free(vrna_fold_compound_t  *vc,
                        vrna_mx_pf_aux_ml_t   *aux_mx)
{
  if (vc && aux_mx) {
    int u;

    free(aux_mx->qqm);
    free(aux_mx->qqm1);

    if (aux_mx->qqmu) {
      for (u = 0; u <= aux_mx->qqmu_size; u++)
        free(aux_mx->qqmu[u]);

      free(aux_mx->qqmu);
    }

    free(aux_mx);
  }
}


PUBLIC FLT_OR_DBL
vrna_exp_E_ml_fast(vrna_fold_compound_t *vc,
                   int                  i,
                   int                  j,
                   vrna_mx_pf_aux_ml_t  *aux_mx)
{
  if (vc) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          return exp_E_ml_fast_window(vc, i, j, aux_mx);
        else
          return exp_E_ml_fast(vc, i, j, aux_mx);

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        return exp_E_ml_fast_comparative(vc, i, j, aux_mx);
        break;

      default:
        vrna_message_warning("vrna_exp_E_ml_fast@multibranch_loops.c: Unknown fold_compound type");
        return 0.;
        break;
    }
  } else {
    return 0.;
  }
}


PRIVATE FLT_OR_DBL
exp_E_ml_fast(vrna_fold_compound_t  *vc,
              int                   i,
              int                   j,
              vrna_mx_pf_aux_ml_t   *aux_mx)
{
  short                     *S1, *S2;
  int                       n, *iidx, k, ij, kl, maxk, ii, with_ud, u, circular, with_gquad,
                            *hc_up_ml, type;
  FLT_OR_DBL                qbt1, temp, *qm, *qb, *qqm, *qqm1, **qqmu, q_temp, q_temp2, *G,
                            *expMLbase,
                            expMLstem;
  vrna_md_t                 *md;
  vrna_exp_param_t          *pf_params;
  vrna_ud_t                 *domains_up;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n                   = (int)vc->length;
  iidx                = vc->iindx;
  ij                  = iidx[i] - j;
  qqm                 = aux_mx->qqm;
  qqm1                = aux_mx->qqm1;
  qqmu                = aux_mx->qqmu;
  qm                  = vc->exp_matrices->qm;
  qb                  = vc->exp_matrices->qb;
  G                   = vc->exp_matrices->G;
  expMLbase           = vc->exp_matrices->expMLbase;
  pf_params           = vc->exp_params;
  md                  = &(pf_params->model_details);
  hc                  = vc->hc;
  sc                  = vc->sc;
  domains_up          = vc->domains_up;
  circular            = md->circ;
  with_gquad          = md->gquad;
  with_ud             = (domains_up && domains_up->exp_energy_cb);
  hc_up_ml            = hc->up_ml;
  hc_dat_local.idx    = vc->jindx;
  hc_dat_local.mx     = hc->matrix;
  hc_dat_local.hc_up  = hc->up_ml;
  hc_dat_local.cp     = vc->cutpoint;

  if (hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default;
  }

  qbt1    = 0;
  q_temp  = 0.;

  qqm[i] = 0.;

  if (with_ud)
    qqmu[0][i] = 0.;

  if (with_gquad)
    expMLstem = exp_E_MLstem(0, -1, -1, pf_params);

  if (evaluate(i, j, i, j - 1, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
    q_temp = qqm1[i] *
             expMLbase[1];

    if (sc) {
      if (sc->exp_energy_up)
        q_temp *= sc->exp_energy_up[j][1];

      if (sc->exp_f)
        q_temp *= sc->exp_f(i, j, i, j - 1, VRNA_DECOMP_ML_ML, sc->data);
    }

    if (with_ud) {
      int cnt;
      for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
        u = domains_up->uniq_motif_size[cnt];
        if (j - u >= i) {
          if (evaluate(i, j, i, j - u, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
            q_temp2 = qqmu[u][i] *
                      domains_up->exp_energy_cb(vc,
                                                j - u + 1,
                                                j,
                                                VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                                domains_up->data) *
                      expMLbase[u];

            if (sc) {
              if (sc->exp_energy_up)
                q_temp2 *= sc->exp_energy_up[j - u + 1][u];

              if (sc->exp_f)
                q_temp2 *= sc->exp_f(i, j, i, j - 1, VRNA_DECOMP_ML_ML, sc->data);
            }

            q_temp += q_temp2;
          }
        }
      }
      qqmu[0][i] += q_temp;
    }

    qqm[i] += q_temp;
  }

  if (evaluate(i, j, i, j, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
    S1    = vc->sequence_encoding;
    S2    = vc->sequence_encoding2;
    type  = get_pair_type_md(S2[i], S2[j], md);

    qbt1 = qb[ij] * exp_E_MLstem(type,
                                 ((i > 1) || circular) ? S1[i - 1] : -1,
                                 ((j < n) || circular) ? S1[j + 1] : -1,
                                 pf_params);
    if (sc)
      if (sc->exp_f)
        qbt1 *= sc->exp_f(i, j, i, j, VRNA_DECOMP_ML_STEM, sc->data);

    qqm[i] += qbt1;

    if (with_ud)
      qqmu[0][i] += qbt1;
  }

  if (with_gquad) {
    /*include gquads into qqm*/
    qqm[i] += G[ij] * expMLstem;

    if (with_ud)
      qqmu[0][i] += G[ij] * expMLstem;
  }

  /*
   *  construction of qm matrix containing multiple loop
   *  partition function contributions from segment i,j
   */
  temp  = 0.0;
  kl    = iidx[i] - j + 1; /* ii-k=[i,k-1] */
  if (hc->f) {
    if (sc && sc->exp_f) {
      for (k = j; k > i; k--, kl++) {
        if (hc->f(i, j, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc->data)) {
          q_temp  = qm[kl] * qqm[k];
          q_temp  *= sc->exp_f(i, j, k - 1, k, VRNA_DECOMP_ML_ML_ML, sc->data);
          temp    += q_temp;
        }
      }
    } else {
      for (k = j; k > i; k--, kl++)
        if (hc->f(i, j, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc->data))
          temp += qm[kl] * qqm[k];
    }
  } else {
    if (sc && sc->exp_f) {
      for (k = j; k > i; k--, kl++) {
        q_temp  = qm[kl] * qqm[k];
        q_temp  *= sc->exp_f(i, j, k - 1, k, VRNA_DECOMP_ML_ML_ML, sc->data);
        temp    += q_temp;
      }
    } else {
      for (k = j; k > i; k--, kl++)
        temp += qm[kl] * qqm[k];
    }
  }

  maxk  = MIN2(i + hc_up_ml[i], j);
  ii    = maxk - i; /* length of unpaired stretch */
  if (with_ud) {
    if (hc->f) {
      if (sc) {
        for (k = maxk; k > i; k--, ii--) {
          if (evaluate(i, j, k, j, VRNA_DECOMP_ML_ML, hc->data)) {
            q_temp = expMLbase[ii] *
                     qqm[k];

            if (sc->exp_energy_up)
              q_temp *= sc->exp_energy_up[i][ii];

            if (sc->exp_f)
              q_temp *= sc->exp_f(i, j, k, j, VRNA_DECOMP_ML_ML, sc->data);

            temp  += q_temp;
            temp  += q_temp *
                     domains_up->exp_energy_cb(vc,
                                               i, k - 1,
                                               VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP,
                                               domains_up->data);
          }
        }
      } else {
        for (k = maxk; k > i; k--, ii--) {
          if (evaluate(i, j, k, j, VRNA_DECOMP_ML_ML, hc->data)) {
            q_temp = expMLbase[ii] *
                     qqm[k];

            temp  += q_temp;
            temp  += q_temp *
                     domains_up->exp_energy_cb(vc,
                                               i, k - 1,
                                               VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP,
                                               domains_up->data);
          }
        }
      }
    } else {
      if (sc) {
        for (k = maxk; k > i; k--, ii--) {
          q_temp = expMLbase[ii] *
                   qqm[k];

          if (sc->exp_energy_up)
            q_temp *= sc->exp_energy_up[i][ii];

          if (sc->exp_f)
            q_temp *= sc->exp_f(i, j, k, j, VRNA_DECOMP_ML_ML, sc->data);

          temp  += q_temp;
          temp  += q_temp *
                   domains_up->exp_energy_cb(vc,
                                             i, k - 1,
                                             VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP,
                                             domains_up->data);
        }
      } else {
        for (k = maxk; k > i; k--, ii--) {
          q_temp = expMLbase[ii] *
                   qqm[k];

          temp  += q_temp;
          temp  += q_temp *
                   domains_up->exp_energy_cb(vc,
                                             i, k - 1,
                                             VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP,
                                             domains_up->data);
        }
      }
    }
  } else {
    if (hc->f) {
      if (sc) {
        for (k = maxk; k > i; k--, ii--) {
          if (hc->f(i, j, k, j, VRNA_DECOMP_ML_ML, hc->data)) {
            q_temp = expMLbase[ii] *
                     qqm[k];

            if (sc->exp_energy_up)
              q_temp *= sc->exp_energy_up[i][ii];

            if (sc->exp_f)
              q_temp *= sc->exp_f(i, j, k, j, VRNA_DECOMP_ML_ML, sc->data);

            temp += q_temp;
          }
        }
      } else {
        for (k = maxk; k > i; k--, ii--)
          if (hc->f(i, j, k, j, VRNA_DECOMP_ML_ML, hc->data))
            temp += expMLbase[ii] *
                    qqm[k];
      }
    } else {
      if (sc) {
        for (k = maxk; k > i; k--, ii--) {
          q_temp = expMLbase[ii] *
                   qqm[k];

          if (sc->exp_energy_up)
            q_temp *= sc->exp_energy_up[i][ii];

          if (sc->exp_f)
            q_temp *= sc->exp_f(i, j, k, j, VRNA_DECOMP_ML_ML, sc->data);

          temp += q_temp;
        }
      } else {
        for (k = maxk; k > i; k--, ii--)
          temp += expMLbase[ii] *
                  qqm[k];
      }
    }
  }

  return temp + qqm[i];
}


PRIVATE FLT_OR_DBL
exp_E_ml_fast_window(vrna_fold_compound_t *vc,
                     int                  i,
                     int                  j,
                     vrna_mx_pf_aux_ml_t  *aux_mx)
{
  short                     *S1;
  int                       n, k, maxk, ii, with_ud, u, circular, with_gquad,
                            *hc_up_ml, type;
  FLT_OR_DBL                qbt1, temp, **qm, **qb, *qqm, *qqm1, **qqmu, q_temp, q_temp2, **G,
                            *expMLbase,
                            expMLstem;
  vrna_md_t                 *md;
  vrna_exp_param_t          *pf_params;
  vrna_ud_t                 *domains_up;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n           = (int)vc->length;
  qqm         = aux_mx->qqm;
  qqm1        = aux_mx->qqm1;
  qqmu        = aux_mx->qqmu;
  qm          = vc->exp_matrices->qm_local;
  qb          = vc->exp_matrices->qb_local;
  G           = vc->exp_matrices->G_local;
  expMLbase   = vc->exp_matrices->expMLbase;
  pf_params   = vc->exp_params;
  md          = &(pf_params->model_details);
  hc          = vc->hc;
  sc          = vc->sc;
  domains_up  = vc->domains_up;
  circular    = md->circ;
  with_gquad  = md->gquad;
  with_ud     = (domains_up && domains_up->exp_energy_cb);
  hc_up_ml    = hc->up_ml;


  hc_dat_local.idx        = vc->jindx;
  hc_dat_local.mx_window  = hc->matrix_local;
  hc_dat_local.hc_up      = hc->up_ml;
  hc_dat_local.cp         = vc->cutpoint;

  if (hc->f) {
    evaluate            = &hc_default_user_window;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default_window;
  }

  qbt1    = 0;
  q_temp  = 0.;

  qqm[i] = 0.;

  if (with_ud)
    qqmu[0][i] = 0.;

  if (with_gquad)
    expMLstem = exp_E_MLstem(0, -1, -1, pf_params);

  if (evaluate(i, j, i, j - 1, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
    q_temp = qqm1[i] *
             expMLbase[1];

    if (sc) {
      if (sc->exp_energy_up)
        q_temp *= sc->exp_energy_up[j][1];

      if (sc->exp_f)
        q_temp *= sc->exp_f(i, j, i, j - 1, VRNA_DECOMP_ML_ML, sc->data);
    }

    if (with_ud) {
      int cnt;
      for (cnt = 0; cnt < domains_up->uniq_motif_count; cnt++) {
        u = domains_up->uniq_motif_size[cnt];
        if (j - u >= i) {
          if (evaluate(i, j, i, j - u, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
            q_temp2 = qqmu[u][i] *
                      domains_up->exp_energy_cb(vc,
                                                j - u + 1,
                                                j,
                                                VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP | VRNA_UNSTRUCTURED_DOMAIN_MOTIF,
                                                domains_up->data) *
                      expMLbase[u];

            if (sc) {
              if (sc->exp_energy_up)
                q_temp2 *= sc->exp_energy_up[j - u + 1][u];

              if (sc->exp_f)
                q_temp2 *= sc->exp_f(i, j, i, j - 1, VRNA_DECOMP_ML_ML, sc->data);
            }

            q_temp += q_temp2;
          }
        }
      }
      qqmu[0][i] += q_temp;
    }

    qqm[i] += q_temp;
  }

  if (evaluate(i, j, i, j, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
    S1    = vc->sequence_encoding;
    type  = get_pair_type_md(S1[i], S1[j], md);

    qbt1 = qb[i][j] * exp_E_MLstem(type,
                                   ((i > 1) || circular) ? S1[i - 1] : -1,
                                   ((j < n) || circular) ? S1[j + 1] : -1,
                                   pf_params);
    if (sc)
      if (sc->exp_f)
        qbt1 *= sc->exp_f(i, j, i, j, VRNA_DECOMP_ML_STEM, sc->data);

    qqm[i] += qbt1;

    if (with_ud)
      qqmu[0][i] += qbt1;
  }

  if (with_gquad) {
    /*include gquads into qqm*/
    qqm[i] += G[i][j] * expMLstem;

    if (with_ud)
      qqmu[0][i] += G[i][j] * expMLstem;
  }

  /*
   *  construction of qm matrix containing multiple loop
   *  partition function contributions from segment i,j
   */
  temp = 0.0;
  if (hc->f) {
    if (sc && sc->exp_f) {
      for (k = j; k > i; k--) {
        if (hc->f(i, j, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc->data)) {
          q_temp  = qm[i][k - 1] * qqm[k];
          q_temp  *= sc->exp_f(i, j, k - 1, k, VRNA_DECOMP_ML_ML_ML, sc->data);
          temp    += q_temp;
        }
      }
    } else {
      for (k = j; k > i; k--)
        if (hc->f(i, j, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc->data))
          temp += qm[i][k - 1] * qqm[k];
    }
  } else {
    if (sc && sc->exp_f) {
      for (k = j; k > i; k--) {
        q_temp  = qm[i][k - 1] * qqm[k];
        q_temp  *= sc->exp_f(i, j, k - 1, k, VRNA_DECOMP_ML_ML_ML, sc->data);
        temp    += q_temp;
      }
    } else {
      for (k = j; k > i; k--)
        temp += qm[i][k - 1] * qqm[k];
    }
  }

  maxk  = MIN2(i + hc_up_ml[i], j);
  ii    = maxk - i; /* length of unpaired stretch */
  if (with_ud) {
    if (hc->f) {
      if (sc) {
        for (k = maxk; k > i; k--, ii--) {
          if (evaluate(i, j, k, j, VRNA_DECOMP_ML_ML, hc->data)) {
            q_temp = expMLbase[ii] *
                     qqm[k];

            if (sc->exp_energy_up)
              q_temp *= sc->exp_energy_up[i][ii];

            if (sc->exp_f)
              q_temp *= sc->exp_f(i, j, k, j, VRNA_DECOMP_ML_ML, sc->data);

            temp  += q_temp;
            temp  += q_temp *
                     domains_up->exp_energy_cb(vc,
                                               i, k - 1,
                                               VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP,
                                               domains_up->data);
          }
        }
      } else {
        for (k = maxk; k > i; k--, ii--) {
          if (evaluate(i, j, k, j, VRNA_DECOMP_ML_ML, hc->data)) {
            q_temp = expMLbase[ii] *
                     qqm[k];

            temp  += q_temp;
            temp  += q_temp *
                     domains_up->exp_energy_cb(vc,
                                               i, k - 1,
                                               VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP,
                                               domains_up->data);
          }
        }
      }
    } else {
      if (sc) {
        for (k = maxk; k > i; k--, ii--) {
          q_temp = expMLbase[ii] *
                   qqm[k];

          if (sc->exp_energy_up)
            q_temp *= sc->exp_energy_up[i][ii];

          if (sc->exp_f)
            q_temp *= sc->exp_f(i, j, k, j, VRNA_DECOMP_ML_ML, sc->data);

          temp  += q_temp;
          temp  += q_temp *
                   domains_up->exp_energy_cb(vc,
                                             i, k - 1,
                                             VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP,
                                             domains_up->data);
        }
      } else {
        for (k = maxk; k > i; k--, ii--) {
          q_temp = expMLbase[ii] *
                   qqm[k];

          temp  += q_temp;
          temp  += q_temp *
                   domains_up->exp_energy_cb(vc,
                                             i, k - 1,
                                             VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP,
                                             domains_up->data);
        }
      }
    }
  } else {
    if (hc->f) {
      if (sc) {
        for (k = maxk; k > i; k--, ii--) {
          if (hc->f(i, j, k, j, VRNA_DECOMP_ML_ML, hc->data)) {
            q_temp = expMLbase[ii] *
                     qqm[k];

            if (sc->exp_energy_up)
              q_temp *= sc->exp_energy_up[i][ii];

            if (sc->exp_f)
              q_temp *= sc->exp_f(i, j, k, j, VRNA_DECOMP_ML_ML, sc->data);

            temp += q_temp;
          }
        }
      } else {
        for (k = maxk; k > i; k--, ii--)
          if (hc->f(i, j, k, j, VRNA_DECOMP_ML_ML, hc->data))
            temp += expMLbase[ii] *
                    qqm[k];
      }
    } else {
      if (sc) {
        for (k = maxk; k > i; k--, ii--) {
          q_temp = expMLbase[ii] *
                   qqm[k];

          if (sc->exp_energy_up)
            q_temp *= sc->exp_energy_up[i][ii];

          if (sc->exp_f)
            q_temp *= sc->exp_f(i, j, k, j, VRNA_DECOMP_ML_ML, sc->data);

          temp += q_temp;
        }
      } else {
        for (k = maxk; k > i; k--, ii--)
          temp += expMLbase[ii] *
                  qqm[k];
      }
    }
  }

  return temp + qqm[i];
}


PRIVATE FLT_OR_DBL
exp_E_ml_fast_comparative(vrna_fold_compound_t  *vc,
                          int                   i,
                          int                   j,
                          vrna_mx_pf_aux_ml_t   *aux_mx)
{
  unsigned int              **a2s;
  short                     **S, **S5, **S3;
  int                       n, s, n_seq, *iidx, k, ij, kl, maxk, ii, circular, *hc_up_ml, type;
  FLT_OR_DBL                temp, *qm, *qb, *qqm, *qqm1, q_temp, *expMLbase;
  vrna_md_t                 *md;
  vrna_exp_param_t          *pf_params;
  vrna_hc_t                 *hc;
  vrna_sc_t                 **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n         = (int)vc->length;
  n_seq     = vc->n_seq;
  iidx      = vc->iindx;
  ij        = iidx[i] - j;
  S         = vc->S;
  S5        = vc->S5;   /* S5[s][i] holds next base 5' of i in sequence s */
  S3        = vc->S3;   /* Sl[s][i] holds next base 3' of i in sequence s */
  a2s       = vc->a2s;
  qqm       = aux_mx->qqm;
  qqm1      = aux_mx->qqm1;
  qm        = vc->exp_matrices->qm;
  qb        = vc->exp_matrices->qb;
  expMLbase = vc->exp_matrices->expMLbase;
  pf_params = vc->exp_params;
  md        = &(pf_params->model_details);
  hc        = vc->hc;
  scs       = vc->scs;
  circular  = md->circ;
  hc_up_ml  = hc->up_ml;

  hc_dat_local.idx    = vc->jindx;
  hc_dat_local.mx     = hc->matrix;
  hc_dat_local.hc_up  = hc->up_ml;
  hc_dat_local.cp     = vc->cutpoint;

  if (hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default;
  }

  q_temp = 0.;

  qqm[i] = 0.;


  if (evaluate(i, j, i, j - 1, VRNA_DECOMP_ML_ML, &hc_dat_local)) {
    q_temp = qqm1[i] *
             expMLbase[1];

    if (scs) {
      for (s = 0; s < n_seq; s++) {
        if (scs[s])
          if (scs[s]->exp_energy_up)
            q_temp *= scs[s]->exp_energy_up[a2s[s][j]][1];
      }
    }

    qqm[i] += q_temp;
  }

  if (evaluate(i, j, i, j, VRNA_DECOMP_ML_STEM, &hc_dat_local)) {
    q_temp = qb[ij];

    for (s = 0; s < n_seq; s++) {
      type    = get_pair_type_md(S[s][i], S[s][j], md);
      q_temp  *= exp_E_MLstem(type,
                              ((i > 1) || circular) ? S5[s][i] : -1,
                              ((j < n) || circular) ? S3[s][j] : -1,
                              pf_params);
    }

    qqm[i] += q_temp;
  }

  /*
   *  construction of qm matrix containing multiple loop
   *  partition function contributions from segment i,j
   */
  temp  = 0.0;
  kl    = iidx[i] - j + 1; /* ii-k=[i,k-1] */
  if (hc->f) {
    for (k = j; k > i; k--, kl++)
      if (hc->f(i, j, k - 1, k, VRNA_DECOMP_ML_ML_ML, hc->data))
        temp += qm[kl] * qqm[k];
  } else {
    for (k = j; k > i; k--, kl++)
      temp += qm[kl] * qqm[k];
  }

  maxk  = MIN2(i + hc_up_ml[i], j);
  ii    = maxk - i; /* length of unpaired stretch */

  if (hc->f) {
    if (scs) {
      for (k = maxk; k > i; k--, ii--) {
        if (hc->f(i, j, k, j, VRNA_DECOMP_ML_ML, hc->data)) {
          q_temp = expMLbase[ii] *
                   qqm[k];

          for (s = 0; s < n_seq; s++) {
            if (scs[s])
              if (scs[s]->exp_energy_up)
                q_temp *= scs[s]->exp_energy_up[a2s[s][i]][a2s[s][k] - a2s[s][i]];
          }
          temp += q_temp;
        }
      }
    } else {
      for (k = maxk; k > i; k--, ii--)
        if (hc->f(i, j, k, j, VRNA_DECOMP_ML_ML, hc->data))
          temp += expMLbase[ii] *
                  qqm[k];
    }
  } else {
    if (scs) {
      for (k = maxk; k > i; k--, ii--) {
        q_temp = expMLbase[ii] *
                 qqm[k];

        for (s = 0; s < n_seq; s++) {
          if (scs[s])
            if (scs[s]->exp_energy_up)
              q_temp *= scs[s]->exp_energy_up[a2s[s][i]][a2s[s][k] - a2s[s][i]];
        }
        temp += q_temp;
      }
    } else {
      for (k = maxk; k > i; k--, ii--)
        temp += expMLbase[ii] *
                qqm[k];
    }
  }

  return temp + qqm[i];
}


PRIVATE unsigned char
hc_default(int            i,
           int            j,
           int            k,
           int            l,
           unsigned char  d,
           void           *data)
{
  int                 ij, kl, di, dj, u;
  unsigned char       eval;
  struct default_data *dat = (struct default_data *)data;

  eval  = (char)0;
  di    = k - i;
  dj    = j - l;

  switch (d) {
    case VRNA_DECOMP_ML_ML_ML:
      u     = l - k - 1;
      eval  = (unsigned char)1;
      if ((u != 0) && (dat->hc_up[k + 1] < u))
        eval = (unsigned char)0;

      break;

    case VRNA_DECOMP_ML_ML:
      eval = (unsigned char)1;
      if ((di != 0) && (dat->hc_up[i] < di))
        eval = (unsigned char)0;

      if ((dj != 0) && (dat->hc_up[l + 1] < dj))
        eval = (unsigned char)0;

      break;

    case VRNA_DECOMP_ML_STEM:
      kl = dat->idx[l] + k;
      if (dat->mx[kl] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) {
        eval = (unsigned char)1;
        if ((di != 0) && (dat->hc_up[i] < di))
          eval = (unsigned char)0;

        if ((dj != 0) && (dat->hc_up[l + 1] < dj))
          eval = (unsigned char)0;
      }

      break;

    case VRNA_DECOMP_PAIR_ML:
      ij = dat->idx[j] + i;
      if (dat->mx[ij] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
        eval = (unsigned char)1;
        di--;
        dj--;
        if ((di != 0) && (dat->hc_up[i + 1] < di))
          eval = (unsigned char)0;

        if ((dj != 0) && (dat->hc_up[l + 1] < dj))
          eval = (unsigned char)0;
      }

      break;

    case VRNA_DECOMP_ML_COAXIAL:
      kl = dat->idx[l] + k;
      if (dat->mx[kl] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC)
        eval = (unsigned char)1;

      break;

    case VRNA_DECOMP_ML_COAXIAL_ENC:
      ij  = dat->idx[j] + i;
      kl  = dat->idx[l] + k;
      if ((dat->mx[ij] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) &&
          (dat->mx[kl] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC))
        eval = (unsigned char)1;

      break;

    default:
      vrna_message_error("hc_cb@multibranch_loops.c: Unrecognized decomposition %d", d);
  }

  return eval;
}


PRIVATE unsigned char
hc_default_window(int           i,
                  int           j,
                  int           k,
                  int           l,
                  unsigned char d,
                  void          *data)
{
  int                 di, dj, u;
  unsigned char       eval;
  struct default_data *dat = (struct default_data *)data;

  eval  = (unsigned char)0;
  di    = k - i;
  dj    = j - l;

  switch (d) {
    case VRNA_DECOMP_ML_ML_ML:
      u     = l - k - 1;
      eval  = (unsigned char)1;
      if ((u != 0) && (dat->hc_up[k + 1] < u))
        eval = (unsigned char)0;

      break;

    case VRNA_DECOMP_ML_ML:
      eval = (unsigned char)1;
      if ((di != 0) && (dat->hc_up[i] < di))
        eval = (unsigned char)0;

      if ((dj != 0) && (dat->hc_up[l + 1] < dj))
        eval = (unsigned char)0;

      break;

    case VRNA_DECOMP_ML_STEM:
      if (dat->mx_window[k][l - k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) {
        eval = (unsigned char)1;
        if ((di != 0) && (dat->hc_up[i] < di))
          eval = (unsigned char)0;

        if ((dj != 0) && (dat->hc_up[l + 1] < dj))
          eval = (unsigned char)0;
      }

      break;

    case VRNA_DECOMP_PAIR_ML:
      if (dat->mx_window[i][j - i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
        eval = (unsigned char)1;
        di--;
        dj--;
        if ((di != 0) && (dat->hc_up[i + 1] < di))
          eval = (unsigned char)0;

        if ((dj != 0) && (dat->hc_up[l + 1] < dj))
          eval = (unsigned char)0;
      }

      break;

    case VRNA_DECOMP_ML_COAXIAL:
      if (dat->mx_window[k][l - k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC)
        eval = (unsigned char)1;

      break;

    case VRNA_DECOMP_ML_COAXIAL_ENC:
      if ((dat->mx_window[i][j - i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC) &&
          (dat->mx_window[k][l - k] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP_ENC))
        eval = (unsigned char)1;

      break;

    default:
      vrna_message_error("hc_cb@multibranch_loops.c: Unrecognized decomposition %d", d);
  }

  return eval;
}


PRIVATE unsigned char
hc_default_ext(int            i,
               int            j,
               int            k,
               int            l,
               unsigned char  d,
               void           *data)
{
  int                 kl, di, dj;
  unsigned char       eval;
  struct default_data *dat = (struct default_data *)data;

  eval  = (unsigned char)0;
  di    = k - i;
  dj    = j - l;
  switch (d) {
    case VRNA_DECOMP_EXT_EXT_STEM:
      kl = dat->idx[j] + l;
      if (dat->mx[kl] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;
        if (i != l) {
          /* otherwise, stem spans from i to j */
          di = l - k - 1;
          if ((di != 0) && (dat->hc_up[k + 1] < di))
            eval = (unsigned char)0;
        }
      }

      break;

    case VRNA_DECOMP_EXT_STEM_EXT:
      kl = dat->idx[k] + i;
      if (dat->mx[kl] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;
        if (j != k) {
          /* otherwise, stem spans from i to j */
          dj = l - k - 1;
          if ((dj != 0) && (dat->hc_up[k + 1] < dj))
            eval = (unsigned char)0;
        }
      }

      break;

    case VRNA_DECOMP_EXT_EXT_STEM1:
      kl = dat->idx[j - 1] + l;
      if (dat->mx[kl] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;

        if (dat->hc_up[j] == 0)
          eval = (unsigned char)0;

        if (i != l) {
          /* otherwise, stem spans from i to j - 1 */
          di = l - k - 1;

          if ((di != 0) && (dat->hc_up[k + 1] < di))
            eval = (unsigned char)0;
        }
      }

      break;

    case VRNA_DECOMP_EXT_STEM_EXT1:
      kl = dat->idx[k] + i + 1;
      if (dat->mx[kl] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;
        if (dat->hc_up[i] == 0)
          eval = (unsigned char)0;

        if (j != k) {
          /* otherwise, stem spans from i + 1 to j */
          dj = l - k - 1;
          if ((dj != 0) && (dat->hc_up[k + 1] < dj))
            eval = (unsigned char)0;
        }
      }

      break;

    case VRNA_DECOMP_EXT_STEM:
      kl = dat->idx[l] + k;
      if (dat->mx[kl] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        eval = (unsigned char)1;
        if ((di != 0) && (dat->hc_up[i] < di))
          eval = (unsigned char)0;

        if ((dj != 0) && (dat->hc_up[l + 1] < dj))
          eval = (unsigned char)0;
      }

      break;

    case VRNA_DECOMP_EXT_EXT:
      eval = (unsigned char)1;
      if ((di != 0) && (dat->hc_up[i] < di))
        eval = (unsigned char)0;

      if ((dj != 0) && (dat->hc_up[l + 1] < dj))
        eval = (unsigned char)0;

      break;

    case VRNA_DECOMP_EXT_UP:
      di    = j - i + 1;
      eval  = (dat->hc_up[i] >= di) ? (unsigned char)1 : (unsigned char)0;
      break;

    default:
      vrna_message_error("hc_cb@multibranch_loops.c: Unrecognized decomposition %d", d);
  }
  return eval;
}


PRIVATE unsigned char
hc_default_user(int           i,
                int           j,
                int           k,
                int           l,
                unsigned char d,
                void          *data)
{
  unsigned char       eval;
  struct default_data *dat = (struct default_data *)data;

  eval  = hc_default(i, j, k, l, d, data);
  eval  = (dat->hc_f(i, j, k, l, d, dat->hc_dat)) ? eval : (unsigned char)0;

  return eval;
}


PRIVATE unsigned char
hc_default_user_window(int            i,
                       int            j,
                       int            k,
                       int            l,
                       unsigned char  d,
                       void           *data)
{
  unsigned char       eval;
  struct default_data *dat = (struct default_data *)data;

  eval  = hc_default_window(i, j, k, l, d, data);
  eval  = (dat->hc_f(i, j, k, l, d, dat->hc_dat)) ? eval : (unsigned char)0;

  return eval;
}


PRIVATE unsigned char
hc_default_user_ext(int           i,
                    int           j,
                    int           k,
                    int           l,
                    unsigned char d,
                    void          *data)
{
  unsigned char       eval;
  struct default_data *dat = (struct default_data *)data;

  eval  = hc_default_ext(i, j, k, l, d, data);
  eval  = (dat->hc_f(i, j, k, l, d, dat->hc_dat)) ? eval : (unsigned char)0;

  return eval;
}
