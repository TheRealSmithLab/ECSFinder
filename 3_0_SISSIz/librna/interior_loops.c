#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include "fold_vars.h"
#include "utils.h"
#include "constraints.h"
#include "exterior_loops.h"
#include "gquad.h"
#include "structured_domains.h"
#include "unstructured_domains.h"
#include "interior_loops.h"


#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

struct default_data {
  void                      *hc_dat;
  vrna_callback_hc_evaluate *hc_f;
};


/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE int
E_int_loop(vrna_fold_compound_t *vc,
           int                  i,
           int                  j);


PRIVATE int
E_int_loop_window(vrna_fold_compound_t  *vc,
                  int                   i,
                  int                   j);


PRIVATE int
E_stack_window(vrna_fold_compound_t *vc,
               int                  i,
               int                  j);


PRIVATE int
E_int_loop_comparative(vrna_fold_compound_t *vc,
                       int                  i,
                       int                  j);


PRIVATE int
E_int_loop_comparative_window(vrna_fold_compound_t  *vc,
                              int                   i,
                              int                   j);


PRIVATE FLT_OR_DBL
exp_E_int_loop(vrna_fold_compound_t *vc,
               int                  i,
               int                  j);


PRIVATE FLT_OR_DBL
exp_E_int_loop_window(vrna_fold_compound_t  *vc,
                      int                   i,
                      int                   j);


PRIVATE FLT_OR_DBL
exp_E_int_loop_comparative(vrna_fold_compound_t *vc,
                           int                  i,
                           int                  j);


PRIVATE INLINE int
eval_interior_loop(vrna_fold_compound_t *vc,
                   int                  i,
                   int                  j,
                   int                  p,
                   int                  q);


PRIVATE INLINE int
eval_int_loop(vrna_fold_compound_t  *vc,
              int                   i,
              int                   j,
              int                   k,
              int                   l);


PRIVATE FLT_OR_DBL
exp_E_interior_loop(vrna_fold_compound_t  *vc,
                    int                   i,
                    int                   j,
                    int                   k,
                    int                   l);


PRIVATE int
BT_int_loop(vrna_fold_compound_t  *vc,
            int                   *i,
            int                   *j,
            int                   en,
            vrna_bp_stack_t       *bp_stack,
            int                   *stack_count);


PRIVATE int
BT_int_loop_comparative(vrna_fold_compound_t  *vc,
                        int                   *i,
                        int                   *j,
                        int                   en,
                        vrna_bp_stack_t       *bp_stack,
                        int                   *stack_count);


PRIVATE int
BT_int_loop_window(vrna_fold_compound_t *vc,
                   int                  *i,
                   int                  *j,
                   int                  en,
                   vrna_bp_stack_t      *bp_stack,
                   int                  *stack_count);


PRIVATE int
BT_int_loop_window_comparative(vrna_fold_compound_t *vc,
                               int                  *i,
                               int                  *j,
                               int                  en,
                               vrna_bp_stack_t      *bp_stack,
                               int                  *stack_count);


PRIVATE int
BT_stack(vrna_fold_compound_t *vc,
         int                  *i,
         int                  *j,
         int                  *en,
         vrna_bp_stack_t      *bp_stack,
         int                  *stack_count);


PRIVATE int
BT_stack_comparative(vrna_fold_compound_t *vc,
                     int                  *i,
                     int                  *j,
                     int                  *en,
                     vrna_bp_stack_t      *bp_stack,
                     int                  *stack_count);


PRIVATE int
BT_stack_window(vrna_fold_compound_t  *vc,
                int                   *i,
                int                   *j,
                int                   *en,
                vrna_bp_stack_t       *bp_stack,
                int                   *stack_count);


PRIVATE int
BT_stack_window_comparative(vrna_fold_compound_t  *vc,
                            int                   *i,
                            int                   *j,
                            int                   *en,
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
hc_default_user(int           i,
                int           j,
                int           k,
                int           l,
                unsigned char d,
                void          *data);


PRIVATE INLINE int get_pair_type(int        i,
                                 int        j,
                                 vrna_md_t  *md);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC int
vrna_E_int_loop(vrna_fold_compound_t  *vc,
                int                   i,
                int                   j)
{
  int e = INF;

  if (vc) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          e = E_int_loop_window(vc, i, j);
        else
          e = E_int_loop(vc, i, j);

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          e = E_int_loop_comparative_window(vc, i, j);
        else
          e = E_int_loop_comparative(vc, i, j);

        break;
    }
  }

  return e;
}


PUBLIC FLT_OR_DBL
vrna_exp_E_int_loop(vrna_fold_compound_t  *vc,
                    int                   i,
                    int                   j)
{
  FLT_OR_DBL q = 0.;

  if (vc) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          q = exp_E_int_loop_window(vc, i, j);
        else
          q = exp_E_int_loop(vc, i, j);

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        q = exp_E_int_loop_comparative(vc, i, j);
        break;
    }
  }

  return q;
}


PUBLIC int
vrna_eval_int_loop(vrna_fold_compound_t *vc,
                   int                  i,
                   int                  j,
                   int                  k,
                   int                  l)
{
  int e = INF;

  if (vc) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        e = eval_int_loop(vc, i, j, k, l);
        break;
    }
  }

  return e;
}


PRIVATE INLINE int
get_pair_type(int       i,
              int       j,
              vrna_md_t *md)
{
  int tt = md->pair[i][j];

  return (tt == 0) ? 7 : tt;
}


PRIVATE INLINE int
eval_int_loop(vrna_fold_compound_t  *vc,
              int                   i,
              int                   j,
              int                   k,
              int                   l)
{
  unsigned int  *sn;
  int           ij, e, cp, *jindx, *rtype, type, type2;
  short         *S, *S2;
  vrna_sc_t     *sc;
  vrna_param_t  *P;
  vrna_md_t     *md;

  cp    = vc->cutpoint;
  jindx = vc->jindx;
  sc    = vc->sc;
  P     = vc->params;
  md    = &(P->model_details);
  sn    = vc->strand_number;
  rtype = &(md->rtype[0]);
  ij    = jindx[j] + i;
  S     = vc->sequence_encoding;
  S2    = vc->sequence_encoding2;

  e = INF;

  if ((sn[k] != sn[i]) || (sn[j] != sn[l]))
    return e;

  type  = get_pair_type(S2[i], S2[j], md);
  type2 = get_pair_type(S2[l], S2[k], md);

  e = ubf_eval_int_loop(i, j, k, l,
                        i + 1, j - 1, k - 1, l + 1,
                        S[i + 1], S[j - 1], S[k - 1], S[l + 1],
                        type, type2, rtype,
                        ij, cp,
                        P, sc);

  return e;
}


PRIVATE int
E_int_loop(vrna_fold_compound_t *vc,
           int                  i,
           int                  j)
{
  unsigned char             type, type_2;
  char                      *ptype, *ptype_pq;
  unsigned char             *hc_pq, *hc, eval_loop;
  short                     *S, S_i1, S_j1, *S_p1, *S_q1;
  unsigned int              *sn;
  int                       q, p, j_q, p_i, pq, *c_pq, max_q, max_p, tmp,
                            *rtype, noGUclosure, no_close, energy, cp,
                            *indx, *hc_up, ij, hc_decompose, e, *c, *ggg,
                            with_gquad, turn;
  vrna_sc_t                 *sc;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_ud_t                 *domains_up;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  cp            = vc->cutpoint;
  indx          = vc->jindx;
  hc            = vc->hc->matrix;
  hc_up         = vc->hc->up_int;
  P             = vc->params;
  ij            = indx[j] + i;
  hc_decompose  = hc[ij];
  e             = INF;
  sn            = vc->strand_number;
  c             = vc->matrices->c;
  ggg           = vc->matrices->ggg;
  md            = &(P->model_details);
  with_gquad    = md->gquad;
  turn          = md->min_loop_size;
  domains_up    = vc->domains_up;

  if (vc->hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = vc->hc->f;
    hc_dat_local.hc_dat = vc->hc->data;
  } else {
    evaluate = &hc_default;
  }

  /* CONSTRAINED INTERIOR LOOP start */
  if (hc_decompose & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
    /* prepare necessary variables */
    rtype       = &(md->rtype[0]);
    noGUclosure = md->noGUclosure;
    max_q       = i + turn + 2;
    max_q       = MAX2(max_q, j - MAXLOOP - 1);

    ptype     = vc->ptype;
    type      = (unsigned char)ptype[ij];
    no_close  = (((type == 3) || (type == 4)) && noGUclosure);
    S         = vc->sequence_encoding;

    S_i1  = S[i + 1];
    S_j1  = S[j - 1];
    sc    = vc->sc;

    if (type == 0)
      type = 7;

    if (domains_up && domains_up->energy_cb) {
      for (q = j - 1; q >= max_q; q--) {
        j_q = j - q - 1;

        if (hc_up[q + 1] < j_q)
          break;

        pq    = indx[q] + i + 1;
        p_i   = 0;
        max_p = i + 1;
        tmp   = i + 1 + MAXLOOP - j_q;
        max_p = MAX2(max_p, tmp);
        tmp   = q - turn;
        max_p = MIN2(max_p, tmp);
        tmp   = i + 1 + hc_up[i + 1];
        max_p = MIN2(max_p, tmp);
        c_pq  = c + pq;

        ptype_pq  = ptype + pq;
        S_p1      = S + i;
        S_q1      = S + q + 1;

        hc_pq = hc + pq;

        for (p = i + 1; p <= max_p; p++) {
          eval_loop = *hc_pq & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC;
          /* discard this configuration if (p,q) is not allowed to be enclosed pair of an interior loop */
          if (eval_loop && evaluate(i, j, p, q, VRNA_DECOMP_PAIR_IL, &hc_dat_local)) {
            energy = *c_pq;
            if (energy != INF) {
              type_2 = rtype[(unsigned char)*ptype_pq];

              if (type_2 == 0)
                type_2 = 7;

              if (noGUclosure)
                if (no_close || (type_2 == 3) || (type_2 == 4))
                  if ((p > i + 1) || (q < j - 1))
                    continue;

              /* continue unless stack */

              energy  += eval_interior_loop(vc, i, j, p, q);
              e       = MIN2(e, energy);
            }
          }

          hc_pq++;    /* get hc[pq + 1] */
          c_pq++;     /* get c[pq + 1] */
          p_i++;      /* increase unpaired region [i+1...p-1] */

          ptype_pq++; /* get ptype[pq + 1] */
          S_p1++;

          pq++;
        } /* end q-loop */
      }   /* end p-loop */
    } else {
      for (q = j - 1; q >= max_q; q--) {
        j_q = j - q - 1;

        if (hc_up[q + 1] < j_q)
          break;

        pq    = indx[q] + i + 1;
        p_i   = 0;
        max_p = i + 1;
        tmp   = i + 1 + MAXLOOP - j_q;
        max_p = MAX2(max_p, tmp);
        tmp   = q - turn;
        max_p = MIN2(max_p, tmp);
        tmp   = i + 1 + hc_up[i + 1];
        max_p = MIN2(max_p, tmp);
        hc_pq = hc + pq;
        c_pq  = c + pq;

        ptype_pq  = ptype + pq;
        S_p1      = S + i;
        S_q1      = S + q + 1;

        for (p = i + 1; p <= max_p; p++) {
          eval_loop = *hc_pq & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC;
          /* discard this configuration if (p,q) is not allowed to be enclosed pair of an interior loop */
          if (eval_loop && evaluate(i, j, p, q, VRNA_DECOMP_PAIR_IL, &hc_dat_local)) {
            energy = *c_pq;
            if (energy != INF) {
              type_2 = rtype[(unsigned char)*ptype_pq];

              if (noGUclosure)
                if (no_close || (type_2 == 3) || (type_2 == 4))
                  if ((p > i + 1) || (q < j - 1))
                    continue;

              /* continue unless stack */

              if (type_2 == 0)
                type_2 = 7;

              energy += ubf_eval_int_loop(i, j, p, q,
                                          i + 1, j - 1, p - 1, q + 1,
                                          S_i1, S_j1, *S_p1, *S_q1,
                                          type, type_2, rtype,
                                          ij, cp,
                                          P, sc);
              e = MIN2(e, energy);
            }
          }

          hc_pq++;    /* get hc[pq + 1] */
          c_pq++;     /* get c[pq + 1] */
          p_i++;      /* increase unpaired region [i+1...p-1] */

          ptype_pq++; /* get ptype[pq + 1] */
          S_p1++;

          pq++;
        } /* end q-loop */
      }   /* end p-loop */
    }

    if (with_gquad) {
      /* include all cases where a g-quadruplex may be enclosed by base pair (i,j) */
      if ((!no_close) && (sn[j] == sn[i])) {
        energy  = E_GQuad_IntLoop(i, j, type, S, ggg, indx, P);
        e       = MIN2(e, energy);
      }
    }
  }

  return e;
}


PRIVATE int
E_int_loop_window(vrna_fold_compound_t  *vc,
                  int                   i,
                  int                   j)
{
  char          **ptype;
  unsigned char hc_decompose, eval_loop;
  short         *S1, si1, sj1, *sp, sp1;
  int           e, p, q, minq, turn, noGUclosure, type, type_2, no_close, energy, *rtype, **c,
                **ggg, with_gquad;

  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_hc_t     *hc;
  vrna_sc_t     *sc;

  e             = INF;
  S1            = vc->sequence_encoding;
  P             = vc->params;
  md            = &(P->model_details);
  c             = vc->matrices->c_local;
  ggg           = vc->matrices->ggg_local;
  hc            = vc->hc;
  sc            = vc->sc;
  ptype         = vc->ptype_local;
  type          = ptype[i][j - i];
  turn          = md->min_loop_size;
  with_gquad    = md->gquad;
  noGUclosure   = md->noGUclosure;
  hc_decompose  = hc->matrix_local[i][j - i];
  rtype         = &(md->rtype[0]);

  if (type == 0)
    type = 7;

  no_close = (((type == 3) || (type == 4)) && noGUclosure);

  if (hc_decompose & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
    int tmp, maxp;
    maxp  = j - 2 - turn;
    tmp   = i + MAXLOOP + 1;
    if (maxp > tmp)
      maxp = tmp;

    tmp = i + 1 + hc->up_int[i + 1];
    if (maxp > tmp)
      maxp = tmp;

    si1 = S1[i + 1];
    sj1 = S1[j - 1];
    sp  = S1 + i;
    for (p = i + 1; p <= maxp; p++, sp++) {
      int u1, u2;
      sp1 = S1[p - 1];
      u1  = p - i - 1;

      minq  = j - i + p - MAXLOOP - 2;
      tmp   = p + 1 + turn;
      if (minq < tmp)
        minq = tmp;

      char          *ptype_p = ptype[p];
      ptype_p -= p;
      int           *c_p = c[p];
      c_p -= p;
      unsigned char *hc_p = hc->matrix_local[p];
      hc_p -= p;

      u2 = j - minq - 1;
      /* seek to minimal q value according to hard constraints */
      while (hc->up_int[minq] < u2) {
        u2--;
        minq++;
        if (minq >= j)
          break;
      }

      /* duplicated code is faster than conditions in loop */
      if (hc->f) {
        for (q = minq; q < j; q++) {
          eval_loop = hc_p[q] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC;
          /* discard this configuration if (p,q) is not allowed to be enclosed pair of an interior loop */
          if (eval_loop && hc->f(i, j, p, q, VRNA_DECOMP_PAIR_IL, hc->data)) {
            energy = c_p[q];

            if (energy < INF) {
              type_2 = ptype_p[q];

              if (type_2 == 0)
                type_2 = 7;

              type_2 = rtype[type_2];

              if ((noGUclosure) && (no_close || (type_2 == 3) || (type_2 == 4)))
                if ((p > i + 1) || (q < j - 1))
                  goto E_int_loop_window_next_q_hc;

              /* continue unless stack */

              energy += E_IntLoop(u1,
                                  u2,
                                  type,
                                  type_2,
                                  si1,
                                  sj1,
                                  sp1,
                                  S1[q + 1],
                                  P);

              /* add soft constraints */
              if (sc) {
                if (sc->energy_up)
                  energy += sc->energy_up[i + 1][u1] +
                            sc->energy_up[q + 1][u2];

                if (sc->energy_bp_local)
                  energy += sc->energy_bp_local[i][j - i];

                if (sc->energy_stack) {
                  if (u1 + u2 == 0) {
                    energy += sc->energy_stack[i] +
                              sc->energy_stack[p] +
                              sc->energy_stack[q] +
                              sc->energy_stack[j];
                  }
                }

                if (sc->f)
                  energy += sc->f(i, j, p, q, VRNA_DECOMP_PAIR_IL, sc->data);
              }

              e = MIN2(e, energy);
            }
          }

E_int_loop_window_next_q_hc:
          u2--;
        } /* end q-loop */
      } else {
        for (q = minq; q < j; q++) {
          /* discard this configuration if (p,q) is not allowed to be enclosed pair of an interior loop */
          if (hc_p[q] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) {
            energy = c_p[q];

            if (energy < INF) {
              type_2 = ptype_p[q];

              if (type_2 == 0)
                type_2 = 7;

              type_2 = rtype[type_2];

              if ((noGUclosure) && (no_close || (type_2 == 3) || (type_2 == 4)))
                if ((p > i + 1) || (q < j - 1))
                  goto E_int_loop_window_next_q;

              /* continue unless stack */

              energy += E_IntLoop(u1,
                                  u2,
                                  type,
                                  type_2,
                                  si1,
                                  sj1,
                                  sp1,
                                  S1[q + 1],
                                  P);

              if (sc) {
                if (sc->energy_up)
                  energy += sc->energy_up[i + 1][u1] +
                            sc->energy_up[q + 1][u2];

                if (sc->energy_bp_local)
                  energy += sc->energy_bp_local[i][j - i];

                if (sc->energy_stack) {
                  if (u1 + u2 == 0) {
                    energy += sc->energy_stack[i] +
                              sc->energy_stack[p] +
                              sc->energy_stack[q] +
                              sc->energy_stack[j];
                  }
                }

                if (sc->f)
                  energy += sc->f(i, j, p, q, VRNA_DECOMP_PAIR_IL, sc->data);
              }

              e = MIN2(e, energy);
            }
          }

E_int_loop_window_next_q:
          u2--;
        } /* end q-loop */
      }   /* end if (hc->f) */
    }     /* end p-loop */
  }

  if (with_gquad) {
    /* include all cases where a g-quadruplex may be enclosed by base pair (i,j) */
    if (!no_close) {
      energy  = E_GQuad_IntLoop_L(i, j, type, S1, ggg, vc->window_size, P);
      e       = MIN2(e, energy);
    }
  }

  return e;
}


PRIVATE INLINE int
ubf_eval_int_loop_comparative(int           col_i,
                              int           col_j,
                              int           col_p,
                              int           col_q,
                              unsigned char type,
                              unsigned char type_2,
                              int           *rtype,
                              int           ij,
                              int           cp,
                              vrna_param_t  *P,
                              short         *SS,
                              short         *S5,
                              short         *S3,
                              unsigned int  *a2s,
                              vrna_sc_t     *sc)
{
  short si, sj, sp, sq;
  int   energy, u1, u2;
  int   i, j, p, q, i1, j1, p1, q1;

  i   = a2s[col_i];
  j   = a2s[col_j];
  p   = a2s[col_p];
  q   = a2s[col_q];
  i1  = a2s[col_i + 1];
  j1  = a2s[col_j - 1];
  p1  = a2s[col_p - 1];
  q1  = a2s[col_q + 1];

  si  = S3[col_i];
  sj  = S5[col_j];
  sp  = S5[col_p];
  sq  = S3[col_q];

  u1  = p1 - i;
  u2  = j1 - q;

  energy = E_IntLoop(u1, u2, type, type_2, si, sj, sp, sq, P);

  /* add soft constraints */
  if (sc) {
    if (sc->energy_up)
      energy += sc->energy_up[i1][u1]
                + sc->energy_up[q1][u2];

    if (sc->energy_bp)
      energy += sc->energy_bp[ij];

    if (sc->energy_stack) {
      if (u1 + u2 == 0) {
        if (SS[col_i] && SS[col_j] && SS[col_p] && SS[col_q]) {
          /* no gap allowed */
          energy += sc->energy_stack[i] +
                    sc->energy_stack[p] +
                    sc->energy_stack[q] +
                    sc->energy_stack[j];
        }
      }
    }

    if (sc->f)
      energy += sc->f(i, j, p, q, VRNA_DECOMP_PAIR_IL, sc->data);
  }

  return energy;
}


PRIVATE int
E_int_loop_comparative(vrna_fold_compound_t *vc,
                       int                  i,
                       int                  j)
{
  unsigned char             type, type_2;
  unsigned char             *hc_pq, *hc, eval_loop;
  unsigned int              **a2s;
  short                     **SS, **S5, **S3, *S_cons;
  int                       q, p, j_q, p_i, u, pq, *c_pq, min_q, max_q, max_p, tmp,
                            *rtype, *types, dangle_model, energy, c0, s, n_seq, cp,
                            *indx, *hc_up, ij, hc_decompose, e, *c, *ggg, with_gquad,
                            turn;
  vrna_sc_t                 *sc, **scs;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  cp            = vc->cutpoint;
  indx          = vc->jindx;
  hc            = vc->hc->matrix;
  hc_up         = vc->hc->up_int;
  P             = vc->params;
  ij            = indx[j] + i;
  hc_decompose  = hc[ij];
  e             = INF;
  c             = vc->matrices->c;
  ggg           = vc->matrices->ggg;
  md            = &(P->model_details);
  with_gquad    = md->gquad;
  turn          = md->min_loop_size;
  dangle_model  = md->dangles;
  types         = NULL;

  if (vc->hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = vc->hc->f;
    hc_dat_local.hc_dat = vc->hc->data;
  } else {
    evaluate = &hc_default;
  }

  /* CONSTRAINED INTERIOR LOOP start */
  if (hc_decompose & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
    SS      = vc->S;
    S5      = vc->S5;     /* S5[s][i] holds next base 5' of i in sequence s */
    S3      = vc->S3;     /* Sl[s][i] holds next base 3' of i in sequence s */
    a2s     = vc->a2s;
    S_cons  = vc->S_cons;
    scs     = vc->scs;
    n_seq   = vc->n_seq;
    types   = (int *)vrna_alloc(sizeof(int) * n_seq);

    for (s = 0; s < n_seq; s++)
      types[s] = get_pair_type(SS[s][i], SS[s][j], md);

    /* prepare necessary variables */
    rtype = &(md->rtype[0]);
    max_q = i + turn + 2;
    max_q = MAX2(max_q, j - MAXLOOP - 1);

    for (q = j - 1; q >= max_q; q--) {
      j_q = j - q - 1;

      if (hc_up[q + 1] < j_q)
        break;

      pq    = indx[q] + i + 1;
      p_i   = 0;
      max_p = i + 1;
      tmp   = i + 1 + MAXLOOP - j_q;
      max_p = MAX2(max_p, tmp);
      tmp   = q - turn;
      max_p = MIN2(max_p, tmp);
      tmp   = i + 1 + hc_up[i + 1];
      max_p = MIN2(max_p, tmp);
      hc_pq = hc + pq;
      c_pq  = c + pq;

      for (p = i + 1; p <= max_p; p++) {
        eval_loop = *hc_pq & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC;
        /* discard this configuration if (p,q) is not allowed to be enclosed pair of an interior loop */
        if (eval_loop && evaluate(i, j, p, q, VRNA_DECOMP_PAIR_IL, &hc_dat_local)) {
          energy = *c_pq;
          if (energy != INF) {
            for (s = 0; s < n_seq; s++) {
              type_2 = get_pair_type(SS[s][q], SS[s][p], md); /* q,p not p,q! */

              sc = (scs && scs[s]) ? scs[s] : NULL;

              energy += ubf_eval_int_loop_comparative(i, j, p, q,
                                                      types[s], type_2, rtype,
                                                      ij, cp,
                                                      P,
                                                      SS[s],
                                                      S5[s],
                                                      S3[s],
                                                      a2s[s],
                                                      sc);
            }
            e = MIN2(e, energy);
          }
        }

        hc_pq++;    /* get hc[pq + 1] */
        c_pq++;     /* get c[pq + 1] */
        p_i++;      /* increase unpaired region [i+1...p-1] */

        pq++;
      } /* end q-loop */
    }   /* end p-loop */


    if (with_gquad) {
      /* include all cases where a g-quadruplex may be enclosed by base pair (i,j) */
      energy = 0;
      for (s = 0; s < n_seq; s++) {
        type = types[s];
        if (dangle_model == 2)
          energy += P->mismatchI[type][S3[s][i]][S5[s][j]];

        if (type > 2)
          energy += P->TerminalAU;
      }
      for (p = i + 2; p < j - VRNA_GQUAD_MIN_BOX_SIZE; p++) {
        u = p - i - 1;
        if (u > MAXLOOP)
          break;

        if (S_cons[p] != 3)
          continue;

        min_q = j - i + p - MAXLOOP - 2;
        c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
        min_q = MAX2(c0, min_q);
        c0    = j - 1;
        max_q = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
        max_q = MIN2(c0, max_q);
        for (q = min_q; q < max_q; q++) {
          if (S_cons[q] != 3)
            continue;

          c0  = energy + ggg[indx[q] + p] + n_seq * P->internal_loop[u + j - q - 1];
          e   = MIN2(e, c0);
        }
      }

      p = i + 1;
      if (S_cons[p] == 3) {
        if (p < j - VRNA_GQUAD_MIN_BOX_SIZE) {
          min_q = j - i + p - MAXLOOP - 2;
          c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
          min_q = MAX2(c0, min_q);
          c0    = j - 3;
          max_q = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
          max_q = MIN2(c0, max_q);
          for (q = min_q; q < max_q; q++) {
            if (S_cons[q] != 3)
              continue;

            c0  = energy + ggg[indx[q] + p] + n_seq * P->internal_loop[j - q - 1];
            e   = MIN2(e, c0);
          }
        }
      }

      q = j - 1;
      if (S_cons[q] == 3) {
        for (p = i + 4; p < j - VRNA_GQUAD_MIN_BOX_SIZE; p++) {
          u = p - i - 1;
          if (u > MAXLOOP)
            break;

          if (S_cons[p] != 3)
            continue;

          c0  = energy + ggg[indx[q] + p] + n_seq * P->internal_loop[u];
          e   = MIN2(e, c0);
        }
      }
    }
  }

  free(types);
  return e;
}


PRIVATE int
E_int_loop_comparative_window(vrna_fold_compound_t  *vc,
                              int                   i,
                              int                   j)
{
  unsigned char type, type_2;
  unsigned char eval_loop;
  unsigned int  **a2s;
  short         **SS, **S5, **S3, *S_cons;
  int           q, p, j_q, p_i, u, min_q, max_q, max_p, tmp,
                *rtype, *types, dangle_model, energy, c0, s, n_seq, cp,
                *hc_up, hc_decompose, e, **c, **ggg, with_gquad,
                turn;
  vrna_sc_t     **scs;
  vrna_param_t  *P;
  vrna_md_t     *md;
  vrna_hc_t     *hc;

  cp            = vc->cutpoint;
  P             = vc->params;
  hc            = vc->hc;
  hc_up         = hc->up_int;
  hc_decompose  = hc->matrix_local[i][j - i];
  e             = INF;
  c             = vc->matrices->c_local;
  ggg           = vc->matrices->ggg_local;
  md            = &(P->model_details);
  with_gquad    = md->gquad;
  turn          = md->min_loop_size;
  dangle_model  = md->dangles;
  types         = NULL;

  /* CONSTRAINED INTERIOR LOOP start */
  if (hc_decompose & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
    SS      = vc->S;
    S5      = vc->S5;     /* S5[s][i] holds next base 5' of i in sequence s */
    S3      = vc->S3;     /* Sl[s][i] holds next base 3' of i in sequence s */
    a2s     = vc->a2s;
    S_cons  = vc->S_cons;
    scs     = vc->scs;
    n_seq   = vc->n_seq;
    types   = (int *)vrna_alloc(sizeof(int) * n_seq);

    for (s = 0; s < n_seq; s++)
      types[s] = get_pair_type(SS[s][i], SS[s][j], md);

    /* prepare necessary variables */
    rtype = &(md->rtype[0]);

    int tmp;
    max_p = j - 2 - turn;
    tmp   = i + MAXLOOP + 1;
    if (max_p > tmp)
      max_p = tmp;

    tmp = i + 1 + hc->up_int[i + 1];
    if (max_p > tmp)
      max_p = tmp;

    if (scs) {
      for (p = i + 1; p <= max_p; p++) {
        min_q = j - i + p - MAXLOOP - 2;
        tmp   = p + 1 + turn;
        if (min_q < tmp)
          min_q = tmp;

        int           *c_p = c[p];
        c_p -= p;
        unsigned char *hc_p = hc->matrix_local[p];
        hc_p -= p;

        j_q = j - min_q - 1;
        /* seek to minimal q value according to hard constraints */
        while (hc->up_int[min_q] < j_q) {
          j_q--;
          min_q++;
          if (min_q >= j)
            break;
        }

        if (hc->f) {
          for (q = min_q; q < j; q++) {
            /* discard this configuration if (p,q) is not allowed to be enclosed pair of an interior loop */
            eval_loop = hc_p[q] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC;
            if (eval_loop && hc->f(i, j, p, q, VRNA_DECOMP_PAIR_IL, hc->data)) {
              energy = c_p[q];

              if (energy < INF) {
                for (s = 0; s < n_seq; s++) {
                  type_2 = get_pair_type(SS[s][q], SS[s][p], md); /* q,p not p,q! */

                  energy += ubf_eval_int_loop_comparative(i, j, p, q,
                                                          types[s], type_2, rtype,
                                                          0, cp,
                                                          P,
                                                          SS[s],
                                                          S5[s],
                                                          S3[s],
                                                          a2s[s],
                                                          scs[s]);
                }
                e = MIN2(e, energy);
              }
            }

            j_q--;
          } /* end q-loop */
        } else {
          for (q = min_q; q < j; q++) {
            /* discard this configuration if (p,q) is not allowed to be enclosed pair of an interior loop */
            if (hc_p[q] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) {
              energy = c_p[q];

              if (energy < INF) {
                for (s = 0; s < n_seq; s++) {
                  type_2 = get_pair_type(SS[s][q], SS[s][p], md); /* q,p not p,q! */

                  energy += ubf_eval_int_loop_comparative(i, j, p, q,
                                                          types[s], type_2, rtype,
                                                          0, cp,
                                                          P,
                                                          SS[s],
                                                          S5[s],
                                                          S3[s],
                                                          a2s[s],
                                                          scs[s]);
                }
                e = MIN2(e, energy);
              }
            }

            j_q--;
          } /* end q-loop */
        }
      }     /* end p-loop */
    } else {
      for (p = i + 1; p <= max_p; p++) {
        min_q = j - i + p - MAXLOOP - 2;
        tmp   = p + 1 + turn;
        if (min_q < tmp)
          min_q = tmp;

        int           *c_p = c[p];
        c_p -= p;
        unsigned char *hc_p = hc->matrix_local[p];
        hc_p -= p;

        j_q = j - min_q - 1;
        /* seek to minimal q value according to hard constraints */
        while (hc->up_int[min_q] < j_q) {
          j_q--;
          min_q++;
          if (min_q >= j)
            break;
        }

        if (hc->f) {
          for (q = min_q; q < j; q++) {
            /* discard this configuration if (p,q) is not allowed to be enclosed pair of an interior loop */
            eval_loop = hc_p[q] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC;
            if (eval_loop && hc->f(i, j, p, q, VRNA_DECOMP_PAIR_IL, hc->data)) {
              energy = c_p[q];

              if (energy < INF) {
                for (s = 0; s < n_seq; s++) {
                  type_2 = get_pair_type(SS[s][q], SS[s][p], md); /* q,p not p,q! */

                  energy += ubf_eval_int_loop_comparative(i, j, p, q,
                                                          types[s], type_2, rtype,
                                                          0, cp,
                                                          P,
                                                          SS[s],
                                                          S5[s],
                                                          S3[s],
                                                          a2s[s],
                                                          NULL);
                }
                e = MIN2(e, energy);
              }
            }

            j_q--;
          } /* end q-loop */
        } else {
          for (q = min_q; q < j; q++) {
            /* discard this configuration if (p,q) is not allowed to be enclosed pair of an interior loop */
            if (hc_p[q] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) {
              energy = c_p[q];

              if (energy < INF) {
                for (s = 0; s < n_seq; s++) {
                  type_2 = get_pair_type(SS[s][q], SS[s][p], md); /* q,p not p,q! */

                  energy += ubf_eval_int_loop_comparative(i, j, p, q,
                                                          types[s], type_2, rtype,
                                                          0, cp,
                                                          P,
                                                          SS[s],
                                                          S5[s],
                                                          S3[s],
                                                          a2s[s],
                                                          NULL);
                }
                e = MIN2(e, energy);
              }
            }

            j_q--;
          } /* end q-loop */
        }
      }     /* end p-loop */
    }

    if (with_gquad) {
      /* include all cases where a g-quadruplex may be enclosed by base pair (i,j) */
      energy = 0;
      for (s = 0; s < n_seq; s++) {
        type = types[s];
        if (dangle_model == 2)
          energy += P->mismatchI[type][S3[s][i]][S5[s][j]];

        if (type > 2)
          energy += P->TerminalAU;
      }
      for (p = i + 2; p < j - VRNA_GQUAD_MIN_BOX_SIZE; p++) {
        u = p - i - 1;
        if (u > MAXLOOP)
          break;

        if (S_cons[p] != 3)
          continue;

        min_q = j - i + p - MAXLOOP - 2;
        c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
        min_q = MAX2(c0, min_q);
        c0    = j - 1;
        max_q = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
        max_q = MIN2(c0, max_q);
        for (q = min_q; q < max_q; q++) {
          if (S_cons[q] != 3)
            continue;

          c0  = energy + ggg[p][q - p] + n_seq * P->internal_loop[u + j - q - 1];
          e   = MIN2(e, c0);
        }
      }

      p = i + 1;
      if (S_cons[p] == 3) {
        if (p < j - VRNA_GQUAD_MIN_BOX_SIZE) {
          min_q = j - i + p - MAXLOOP - 2;
          c0    = p + VRNA_GQUAD_MIN_BOX_SIZE - 1;
          min_q = MAX2(c0, min_q);
          c0    = j - 3;
          max_q = p + VRNA_GQUAD_MAX_BOX_SIZE + 1;
          max_q = MIN2(c0, max_q);
          for (q = min_q; q < max_q; q++) {
            if (S_cons[q] != 3)
              continue;

            c0  = energy + ggg[p][q - p] + n_seq * P->internal_loop[j - q - 1];
            e   = MIN2(e, c0);
          }
        }
      }

      q = j - 1;
      if (S_cons[q] == 3) {
        for (p = i + 4; p < j - VRNA_GQUAD_MIN_BOX_SIZE; p++) {
          u = p - i - 1;
          if (u > MAXLOOP)
            break;

          if (S_cons[p] != 3)
            continue;

          c0  = energy + ggg[p][q - p] + n_seq * P->internal_loop[u];
          e   = MIN2(e, c0);
        }
      }
    }
  }

  free(types);
  return e;
}


PRIVATE INLINE int
eval_interior_loop(vrna_fold_compound_t *vc,
                   int                  i,
                   int                  j,
                   int                  p,
                   int                  q)
{
  int           energy, e, e5, e3, u1, u2, i1, j1, p1, q1, *idx;
  short         *S, *S2, si, sj, sp, sq;
  unsigned char type, type_2;
  unsigned int  *sn;
  int           *rtype, ij, cp;
  vrna_param_t  *P;
  vrna_sc_t     *sc;
  vrna_md_t     *md;
  vrna_ud_t     *domains_up;

  S           = vc->sequence_encoding;
  S2          = vc->sequence_encoding2;
  cp          = vc->cutpoint;
  P           = vc->params;
  sc          = vc->sc;
  domains_up  = vc->domains_up;
  md          = &(P->model_details);
  sn          = vc->strand_number;
  rtype       = &(md->rtype[0]);
  idx         = vc->jindx;
  type        = get_pair_type(S2[i], S2[j], md);
  type_2      = get_pair_type(S2[q], S2[p], md);

  i1  = i + 1;
  j1  = j - 1;
  p1  = p - 1;
  q1  = q + 1;
  u1  = p1 - i;
  u2  = j1 - q;

  si  = S[i1];
  sj  = S[j1];
  sp  = S[p1];
  sq  = S[q1];
  ij  = idx[j] + i;

  if ((sn[p] == sn[i]) && (sn[j] == sn[q])) {
    /* regular interior loop */
    energy = E_IntLoop(u1, u2, type, type_2, si, sj, sp, sq, P);
  } else {
    /* interior loop like cofold structure */
    short Si, Sj;
    Si      = (sn[i1] == sn[i]) ? si : -1;
    Sj      = (sn[j] == sn[j1]) ? sj : -1;
    energy  = E_IntLoop_Co(rtype[type], rtype[type_2],
                           i, j, p, q,
                           cp,
                           Si, Sj,
                           sp, sq,
                           P->model_details.dangles,
                           P);
  }

  /* add soft constraints */
  if (sc) {
    if (sc->energy_up)
      energy += sc->energy_up[i1][u1] +
                sc->energy_up[q1][u2];

    if (sc->energy_bp)
      energy += sc->energy_bp[ij];

    if (sc->energy_stack) {
      if (u1 + u2 == 0) {
        energy += sc->energy_stack[i] +
                  sc->energy_stack[p] +
                  sc->energy_stack[q] +
                  sc->energy_stack[j];
      }
    }

    if (sc->f)
      energy += sc->f(i, j, p, q, VRNA_DECOMP_PAIR_IL, sc->data);
  }

  e   = energy;
  e5  = e3 = 0;

  if (u1 > 0) {
    e5 = domains_up->energy_cb(vc,
                               i + 1, p - 1,
                               VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                               domains_up->data);
  }

  if (u2 > 0) {
    e3 = domains_up->energy_cb(vc,
                               q + 1, j - 1,
                               VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                               domains_up->data);
  }

  e = MIN2(e, energy + e5);
  e = MIN2(e, energy + e3);
  e = MIN2(e, energy + e5 + e3);


  return e;
}


PRIVATE FLT_OR_DBL
exp_E_int_loop(vrna_fold_compound_t *vc,
               int                  i,
               int                  j)
{
  unsigned char             type, type_2;
  char                      *ptype;
  unsigned char             *hc, eval_loop;
  short                     *S1, S_i1, S_j1;
  unsigned int              *sn;
  int                       k, l, u1, u2, kl, maxk, minl, *rtype, noGUclosure,
                            no_close, *my_iindx, *jindx, *hc_up, ij,
                            with_gquad, turn;
  FLT_OR_DBL                qbt1, q_temp, *qb, *G, *scale;
  vrna_sc_t                 *sc;
  vrna_exp_param_t          *pf_params;
  vrna_md_t                 *md;
  vrna_ud_t                 *domains_up;
  vrna_callback_hc_evaluate *evaluate;
  struct  default_data      hc_dat_local;

  ptype       = vc->ptype;
  S1          = vc->sequence_encoding;
  S_i1        = S1[i + 1];
  S_j1        = S1[j - 1];
  my_iindx    = vc->iindx;
  jindx       = vc->jindx;
  hc          = vc->hc->matrix;
  hc_up       = vc->hc->up_int;
  sc          = vc->sc;
  sn          = vc->strand_number;
  pf_params   = vc->exp_params;
  ij          = jindx[j] + i;
  md          = &(pf_params->model_details);
  with_gquad  = md->gquad;
  turn        = md->min_loop_size;
  qb          = vc->exp_matrices->qb;
  G           = vc->exp_matrices->G;
  scale       = vc->exp_matrices->scale;
  domains_up  = vc->domains_up;
  qbt1        = 0.;

  if (vc->hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = vc->hc->f;
    hc_dat_local.hc_dat = vc->hc->data;
  } else {
    evaluate = &hc_default;
  }

  /* CONSTRAINED INTERIOR LOOP start */
  if (hc[ij] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
    type        = (unsigned char)ptype[ij];
    rtype       = &(md->rtype[0]);
    noGUclosure = md->noGUclosure;
    no_close    = (((type == 3) || (type == 4)) && noGUclosure);
    maxk        = i + MAXLOOP + 1;
    maxk        = MIN2(maxk, j - turn - 2);
    maxk        = MIN2(maxk, i + 1 + hc_up[i + 1]);

    if (type == 0)
      type = 7;

    for (k = i + 1; k <= maxk; k++) {
      if (sn[k] != sn[i])
        break;

      u1 = k - i - 1;

      minl  = MAX2(k + turn + 1, j - 1 - MAXLOOP + u1);
      kl    = my_iindx[k] - j + 1;

      for (u2 = 0, l = j - 1; l >= minl; l--, kl++, u2++) {
        if (hc_up[l + 1] < u2)
          break;

        eval_loop =
          (hc[jindx[l] + k] &
           VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) ? (unsigned char)1 : (unsigned char)0;

        /* discard this configuration if (p,q) is not allowed to be enclosed pair of an interior loop */
        if (eval_loop && evaluate(i, j, k, l, VRNA_DECOMP_PAIR_IL, &hc_dat_local)) {
          if (sn[j] != sn[l])
            break;

          type_2 = rtype[(unsigned char)ptype[jindx[l] + k]];

          if (type_2 == 0)
            type_2 = 7;

          q_temp = qb[kl] *
                   scale[u1 + u2 + 2] *
                   exp_E_IntLoop(u1, u2, type, type_2, S_i1, S_j1, S1[k - 1], S1[l + 1], pf_params);

          /* soft constraints */
          if (sc) {
            if (sc->exp_energy_up)
              q_temp *= sc->exp_energy_up[i + 1][u1] *
                        sc->exp_energy_up[l + 1][u2];

            if (sc->exp_f)
              q_temp *= sc->exp_f(i, j, k, l, VRNA_DECOMP_PAIR_IL, sc->data);

            if (sc->exp_energy_stack) {
              if ((i + 1 == k) && (j - 1 == l)) {
                q_temp *= sc->exp_energy_stack[i] *
                          sc->exp_energy_stack[k] *
                          sc->exp_energy_stack[l] *
                          sc->exp_energy_stack[j];
              }
            }
          }

          qbt1 += q_temp;

          /* unstructured domains */
          if (domains_up && domains_up->exp_energy_cb) {
            FLT_OR_DBL qq5, qq3;

            qq5 = qq3 = 0.;

            if (u1 > 0) {
              qq5 = domains_up->exp_energy_cb(vc,
                                              i + 1, k - 1,
                                              VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                              domains_up->data);
            }

            if (u2 > 0) {
              qq3 = domains_up->exp_energy_cb(vc,
                                              l + 1, j - 1,
                                              VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                              domains_up->data);
            }

            qbt1  += q_temp * qq5;        /* only motifs in 5' part */
            qbt1  += q_temp * qq3;        /* only motifs in 3' part */
            qbt1  += q_temp * qq5 * qq3;  /* motifs in both parts */
          }
        }
      }
    }

    if (with_gquad) {
      /* include all cases where a g-quadruplex may be enclosed by base pair (i,j) */
      if ((!no_close) && (sn[j] == sn[i]))
        qbt1 += exp_E_GQuad_IntLoop(i, j, type, S1, G, scale, my_iindx, pf_params);
    }

    if (sc && sc->exp_energy_bp)
      qbt1 *= sc->exp_energy_bp[my_iindx[i] - j];
  }

  return qbt1;
}


PRIVATE FLT_OR_DBL
exp_E_int_loop_window(vrna_fold_compound_t  *vc,
                      int                   i,
                      int                   j)
{
  unsigned char             type, type_2;
  char                      **ptype;
  unsigned char             **hc, eval_loop;
  short                     *S1, S_i1, S_j1;
  unsigned int              *sn;
  int                       k, l, u1, u2, maxk, minl, *rtype, noGUclosure,
                            no_close, *my_iindx, *jindx, *hc_up,
                            with_gquad, turn;
  FLT_OR_DBL                qbt1, q_temp, **qb, **G, *scale;
  vrna_sc_t                 *sc;
  vrna_exp_param_t          *pf_params;
  vrna_md_t                 *md;
  vrna_ud_t                 *domains_up;
  vrna_callback_hc_evaluate *evaluate;
  struct  default_data      hc_dat_local;

  ptype       = vc->ptype_local;
  S1          = vc->sequence_encoding;
  S_i1        = S1[i + 1];
  S_j1        = S1[j - 1];
  jindx       = vc->jindx;
  hc          = vc->hc->matrix_local;
  hc_up       = vc->hc->up_int;
  sc          = vc->sc;
  sn          = vc->strand_number;
  pf_params   = vc->exp_params;
  md          = &(pf_params->model_details);
  with_gquad  = md->gquad;
  turn        = md->min_loop_size;
  qb          = vc->exp_matrices->qb_local;
  G           = vc->exp_matrices->G_local;
  scale       = vc->exp_matrices->scale;
  domains_up  = vc->domains_up;
  qbt1        = 0.;

  if (vc->hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = vc->hc->f;
    hc_dat_local.hc_dat = vc->hc->data;
  } else {
    evaluate = &hc_default;
  }

  /* CONSTRAINED INTERIOR LOOP start */
  if (hc[i][j - i] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
    type        = ptype[i][j];
    rtype       = &(md->rtype[0]);
    noGUclosure = md->noGUclosure;
    no_close    = (((type == 3) || (type == 4)) && noGUclosure);
    maxk        = i + MAXLOOP + 1;
    maxk        = MIN2(maxk, j - turn - 2);
    maxk        = MIN2(maxk, i + 1 + hc_up[i + 1]);

    if (type == 0)
      type = 7;

    for (k = i + 1; k <= maxk; k++) {
      if (sn[k] != sn[i])
        break;

      u1 = k - i - 1;

      minl = MAX2(k + turn + 1, j - 1 - MAXLOOP + u1);

      for (u2 = 0, l = j - 1; l >= minl; l--, u2++) {
        if (hc_up[l + 1] < u2)
          break;

        eval_loop =
          (hc[k][l - k] &
           VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) ? (unsigned char)1 : (unsigned char)0;

        /* discard this configuration if (p,q) is not allowed to be enclosed pair of an interior loop */
        if (eval_loop && evaluate(i, j, k, l, VRNA_DECOMP_PAIR_IL, &hc_dat_local)) {
          if (sn[j] != sn[l])
            break;

          type_2 = rtype[(unsigned char)ptype[k][l]];

          if (type_2 == 0)
            type_2 = 7;

          q_temp = qb[k][l] *
                   scale[u1 + u2 + 2] *
                   exp_E_IntLoop(u1, u2, type, type_2, S_i1, S_j1, S1[k - 1], S1[l + 1], pf_params);

          /* soft constraints */
          if (sc) {
            if (sc->exp_energy_up)
              q_temp *= sc->exp_energy_up[i + 1][u1] *
                        sc->exp_energy_up[l + 1][u2];

            if (sc->exp_energy_bp_local)
              q_temp *= sc->exp_energy_bp_local[i][j - i];

            if (sc->exp_f)
              q_temp *= sc->exp_f(i, j, k, l, VRNA_DECOMP_PAIR_IL, sc->data);

            if (sc->exp_energy_stack) {
              if ((u1 == 0) && (u2 == 0)) {
                q_temp *= sc->exp_energy_stack[i] *
                          sc->exp_energy_stack[k] *
                          sc->exp_energy_stack[l] *
                          sc->exp_energy_stack[j];
              }
            }
          }

          qbt1 += q_temp;

          /* unstructured domains */
          if (domains_up && domains_up->exp_energy_cb) {
            FLT_OR_DBL qq5, qq3;

            qq5 = qq3 = 0.;

            if (u1 > 0) {
              qq5 = domains_up->exp_energy_cb(vc,
                                              i + 1, k - 1,
                                              VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                              domains_up->data);
            }

            if (u2 > 0) {
              qq3 = domains_up->exp_energy_cb(vc,
                                              l + 1, j - 1,
                                              VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                              domains_up->data);
            }

            qbt1  += q_temp * qq5;        /* only motifs in 5' part */
            qbt1  += q_temp * qq3;        /* only motifs in 3' part */
            qbt1  += q_temp * qq5 * qq3;  /* motifs in both parts */
          }
        }
      }
    }

#if 0
    /* no G-Quadruplexes for sliding-window partition function yet! */
    if (with_gquad) {
      /* include all cases where a g-quadruplex may be enclosed by base pair (i,j) */
      if ((!no_close) && (sn[j] == sn[i]))
        qbt1 += exp_E_GQuad_IntLoop(i, j, type, S1, G, scale, my_iindx, pf_params);
    }

#endif
  }

  return qbt1;
}


PUBLIC FLT_OR_DBL
vrna_exp_E_interior_loop(vrna_fold_compound_t *vc,
                         int                  i,
                         int                  j,
                         int                  k,
                         int                  l)
{
  if (vc)
    return exp_E_interior_loop(vc, i, j, k, l);

  return 0.;
}


PRIVATE FLT_OR_DBL
exp_E_interior_loop(vrna_fold_compound_t  *vc,
                    int                   i,
                    int                   j,
                    int                   k,
                    int                   l)
{
  unsigned char             type, type_2;
  char                      *ptype;
  unsigned char             *hc, eval_loop;
  short                     *S1, S_i1, S_j1;
  unsigned int              *sn;
  int                       u1, u2, *rtype, *my_iindx, *jindx, *hc_up, ij;
  FLT_OR_DBL                qbt1, q_temp, *scale;
  vrna_sc_t                 *sc;
  vrna_exp_param_t          *pf_params;
  vrna_md_t                 *md;
  vrna_ud_t                 *domains_up;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  ptype       = vc->ptype;
  S1          = vc->sequence_encoding;
  S_i1        = S1[i + 1];
  S_j1        = S1[j - 1];
  my_iindx    = vc->iindx;
  jindx       = vc->jindx;
  hc          = vc->hc->matrix;
  hc_up       = vc->hc->up_int;
  sc          = vc->sc;
  pf_params   = vc->exp_params;
  ij          = jindx[j] + i;
  sn          = vc->strand_number;
  md          = &(pf_params->model_details);
  scale       = vc->exp_matrices->scale;
  domains_up  = vc->domains_up;
  qbt1        = 0.;
  u1          = k - i - 1;
  u2          = j - l - 1;

  if ((sn[k] != sn[i]) || (sn[j] != sn[l]))
    return qbt1;

  if (hc_up[l + 1] < u2)
    return qbt1;

  if (hc_up[i + 1] < u1)
    return qbt1;

  if (vc->hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = vc->hc->f;
    hc_dat_local.hc_dat = vc->hc->data;
  } else {
    evaluate = &hc_default;
  }

  /* CONSTRAINED INTERIOR LOOP start */
  eval_loop =
    ((hc[ij] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) &&
     (hc[jindx[l] + k] &
      VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC)) ? (unsigned char)1 : (unsigned char)0;

  /* discard this configuration if (p,q) is not allowed to be enclosed pair of an interior loop */
  if (eval_loop && evaluate(i, j, k, l, VRNA_DECOMP_PAIR_IL, &hc_dat_local)) {
    type    = (unsigned char)ptype[ij];
    rtype   = &(md->rtype[0]);
    type    = (unsigned char)ptype[ij];
    type_2  = rtype[(unsigned char)ptype[jindx[l] + k]];

    if (type == 0)
      type = 7;

    if (type_2 == 0)
      type_2 = 7;

    q_temp = exp_E_IntLoop(u1, u2, type, type_2, S_i1, S_j1, S1[k - 1], S1[l + 1], pf_params) *
             scale[u1 + u2 + 2];

    /* soft constraints */
    if (sc) {
      if (sc->exp_energy_up)
        q_temp *= sc->exp_energy_up[i + 1][u1] *
                  sc->exp_energy_up[l + 1][u2];

      if (sc->exp_f)
        q_temp *= sc->exp_f(i, j, k, l, VRNA_DECOMP_PAIR_IL, sc->data);

      if (sc->exp_energy_stack) {
        if ((i + 1 == k) && (j - 1 == l)) {
          q_temp *= sc->exp_energy_stack[i] *
                    sc->exp_energy_stack[k] *
                    sc->exp_energy_stack[l] *
                    sc->exp_energy_stack[j];
        }
      }

      if (sc->exp_energy_bp)
        q_temp *= sc->exp_energy_bp[my_iindx[i] - j];
    }

    qbt1 += q_temp;

    /* unstructured domains */
    if (domains_up && domains_up->exp_energy_cb) {
      FLT_OR_DBL qq5, qq3;

      qq5 = qq3 = 0.;

      if (u1 > 0) {
        qq5 = domains_up->exp_energy_cb(vc,
                                        i + 1, k - 1,
                                        VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                        domains_up->data);
      }

      if (u2 > 0) {
        qq3 = domains_up->exp_energy_cb(vc,
                                        l + 1, j - 1,
                                        VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                                        domains_up->data);
      }

      qbt1  += q_temp * qq5;        /* only motifs in 5' part */
      qbt1  += q_temp * qq3;        /* only motigs in 3' part */
      qbt1  += q_temp * qq5 * qq3;  /* motifs in both parts */
    }
  }

  return qbt1;
}


PRIVATE FLT_OR_DBL
exp_E_int_loop_comparative(vrna_fold_compound_t *vc,
                           int                  i,
                           int                  j)
{
  unsigned char             type_2;
  unsigned char             *hc, eval_loop;
  unsigned int              **a2s;
  short                     **S, **S5, **S3;
  int                       n_seq, s, ij, jij, k, l, u1, u2, kl, maxk, minl, *types,
                            turn, with_gquad, *hc_up, *jindx, *my_iindx;
  FLT_OR_DBL                qbt1, *qb, *scale, qloop;
  vrna_sc_t                 **scs;
  vrna_exp_param_t          *pf_params;
  vrna_md_t                 *md;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  types       = NULL;
  my_iindx    = vc->iindx;
  jindx       = vc->jindx;
  hc          = vc->hc->matrix;
  hc_up       = vc->hc->up_int;
  pf_params   = vc->exp_params;
  md          = &(pf_params->model_details);
  with_gquad  = md->gquad;
  turn        = md->min_loop_size;
  qb          = vc->exp_matrices->qb;
  scale       = vc->exp_matrices->scale;
  qbt1        = 0.;
  jij         = jindx[j] + i;
  ij          = my_iindx[i] - j;

  if (vc->hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = vc->hc->f;
    hc_dat_local.hc_dat = vc->hc->data;
  } else {
    evaluate = &hc_default;
  }

  /* CONSTRAINED INTERIOR LOOP start */
  if (hc[jij] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
    S     = vc->S;
    S5    = vc->S5;     /* S5[s][i] holds next base 5' of i in sequence s */
    S3    = vc->S3;     /* Sl[s][i] holds next base 3' of i in sequence s */
    a2s   = vc->a2s;
    scs   = vc->scs;
    n_seq = vc->n_seq;
    types = (int *)vrna_alloc(sizeof(int) * n_seq);

    for (s = 0; s < n_seq; s++)
      types[s] = get_pair_type(S[s][i], S[s][j], md);

    /* prepare necessary variables */
    maxk  = i + MAXLOOP + 1;
    maxk  = MIN2(maxk, j - turn - 2);
    maxk  = MIN2(maxk, i + 1 + hc_up[i + 1]);

    for (k = i + 1; k <= maxk; k++) {
      u1 = k - i - 1;

      minl  = MAX2(k + turn + 1, j - 1 - MAXLOOP + u1);
      kl    = my_iindx[k] - j + 1;

      for (l = j - 1; l >= minl; l--, kl++, u2++) {
        if (hc_up[l + 1] < j - l - 1)
          break;

        eval_loop =
          (hc[jindx[l] + k] &
           VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) ? (unsigned char)1 : (unsigned char)0;

        /* discard this configuration if (p,q) is not allowed to be enclosed pair of an interior loop */
        if (eval_loop && evaluate(i, j, k, l, VRNA_DECOMP_PAIR_IL, &hc_dat_local)) {
          qloop = 1.;

          for (s = 0; s < n_seq; s++) {
            u1      = a2s[s][k - 1] - a2s[s][i];
            u2      = a2s[s][j - 1] - a2s[s][l];
            type_2  = get_pair_type(S[s][l], S[s][k], md);

            qloop *= exp_E_IntLoop(u1, u2,
                                   types[s], type_2, S3[s][i],
                                   S5[s][j], S5[s][k], S3[s][l],
                                   pf_params
                                   );
          }

          if (scs) {
            for (s = 0; s < n_seq; s++) {
              if (scs[s]) {
                u1  = a2s[s][k - 1] - a2s[s][i];
                u2  = a2s[s][j - 1] - a2s[s][l];

                if (scs[s]->exp_energy_up)
                  qloop *= scs[s]->exp_energy_up[a2s[s][i] + 1][u1] *
                           scs[s]->exp_energy_up[a2s[s][l] + 1][u2];

                if (scs[s]->exp_energy_stack) {
                  if (u1 + u2 == 0) {
                    if (S[s][i] && S[s][j] && S[s][k] && S[s][l]) {
                      /* don't allow gaps in stack */
                      qloop *= scs[s]->exp_energy_stack[i] *
                               scs[s]->exp_energy_stack[k] *
                               scs[s]->exp_energy_stack[l] *
                               scs[s]->exp_energy_stack[j];
                    }
                  }
                }
              }
            }
          }

          qbt1 += qb[my_iindx[k] - l] *
                  qloop *
                  scale[k - i + j - l];
        }
      }
    }

    if (with_gquad) {
      /* include all cases where a g-quadruplex may be enclosed by base pair (i,j) */
      /* not implemented yet! */
    }

    if (scs) {
      for (s = 0; s < n_seq; s++)
        if (scs[s] && scs[s]->exp_energy_bp)
          qbt1 *= scs[s]->exp_energy_bp[ij];
    }
  }

  /* cleanup */
  free(types);

  return qbt1;
}


PUBLIC int
vrna_E_ext_int_loop(vrna_fold_compound_t  *vc,
                    int                   i,
                    int                   j,
                    int                   *ip,
                    int                   *iq)
{
  unsigned char             type, type_2;
  int                       ij, q, p, e, s, u1, u2, qmin, energy, *rtype, *types,
                            length, *indx, *hc_up, *c, turn, n_seq;
  char                      *ptype;
  unsigned char             *hc, eval_loop;
  unsigned int              **a2s;
  short                     *S, **SS, **S5, **S3;
  vrna_md_t                 *md;
  vrna_param_t              *P;
  vrna_sc_t                 *sc, **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  length  = vc->length;
  indx    = vc->jindx;
  ptype   = vc->ptype;
  c       = vc->matrices->c;
  hc      = vc->hc->matrix;
  hc_up   = vc->hc->up_int;
  P       = vc->params;
  md      = &(P->model_details);
  turn    = md->min_loop_size;
  types   = NULL;
  ij      = indx[j] + i;
  rtype   = &(md->rtype[0]);
  e       = INF;

  if (vc->hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = vc->hc->f;
    hc_dat_local.hc_dat = vc->hc->data;
  } else {
    evaluate = &hc_default;
  }

  /* CONSTRAINED INTERIOR LOOP start */
  if (hc[ij] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
    /* prepare necessary variables */
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        type = rtype[(unsigned char)ptype[ij]];

        if (type == 0)
          type = 7;

        S   = vc->sequence_encoding;
        sc  = vc->sc;
        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        SS    = vc->S;
        S5    = vc->S5;   /* S5[s][i] holds next base 5' of i in sequence s */
        S3    = vc->S3;   /* Sl[s][i] holds next base 3' of i in sequence s */
        a2s   = vc->a2s;
        scs   = vc->scs;
        n_seq = vc->n_seq;
        types = (int *)vrna_alloc(sizeof(int) * n_seq);

        for (s = 0; s < n_seq; s++)
          types[s] = get_pair_type(SS[s][j], SS[s][i], md);
        break;

      default:
        return e;
        break;
    }

    for (p = j + 1; p < length; p++) {
      u1 = p - j - 1;
      if (u1 + i - 1 > MAXLOOP)
        break;

      if (hc_up[j + 1] < u1)
        break;

      qmin = u1 + i - 1 + length - MAXLOOP;
      if (qmin < p + turn + 1)
        qmin = p + turn + 1;

      for (q = length; q >= qmin; q--) {
        u2 = i - 1 + length - q;
        if (hc_up[q + 1] < u2)
          break;

        if (u1 + u2 > MAXLOOP)
          continue;

        eval_loop = hc[indx[q] + p] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP;

        if (eval_loop && evaluate(i, j, p, q, VRNA_DECOMP_PAIR_IL, &hc_dat_local)) {
          energy = c[indx[q] + p];
          if (energy < INF) {
            switch (vc->type) {
              case VRNA_FC_TYPE_SINGLE:
                type_2 = rtype[(unsigned char)ptype[indx[q] + p]];

                if (type_2 == 0)
                  type_2 = 7;

                energy += ubf_eval_ext_int_loop(i, j, p, q,
                                                i - 1, j + 1, p - 1, q + 1,
                                                S[j + 1], S[i - 1], S[p - 1], S[q + 1],
                                                type, type_2,
                                                length,
                                                P, sc);
                break;

              case VRNA_FC_TYPE_COMPARATIVE:
                energy = 0;
                for (s = 0; s < n_seq; s++) {
                  type_2 = get_pair_type(SS[s][q], SS[s][p], md); /* q,p not p,q! */

                  sc = (scs && scs[s]) ? scs[s] : NULL;

                  energy += ubf_eval_ext_int_loop(a2s[s][i],
                                                  a2s[s][j],
                                                  a2s[s][p],
                                                  a2s[s][q],
                                                  a2s[s][i - 1],
                                                  a2s[s][j + 1],
                                                  a2s[s][p - 1],
                                                  a2s[s][q + 1],
                                                  S3[s][j],
                                                  S5[s][i],
                                                  S5[s][p],
                                                  S3[s][q],
                                                  types[s],
                                                  type_2,
                                                  a2s[s][length],
                                                  P,
                                                  sc);
                }
                break;
            }

            if (energy < e) {
              e = energy;
              if ((ip != NULL) && (iq != NULL)) {
                *ip = p;
                *iq = q;
              }
            }
          }
        }
      }
    }
  }

  free(types);

  return e;
}


PUBLIC int
vrna_E_stack(vrna_fold_compound_t *vc,
             int                  i,
             int                  j)
{
  unsigned char             type, type_2;
  char                      *ptype;
  unsigned char             *hard_constraints, eval_loop;
  unsigned int              **a2s;
  short                     *S, **SS;
  unsigned int              *sn;
  int                       e, ij, pq, p, q, s, n_seq, cp, *rtype, *indx;
  vrna_sc_t                 *sc, **scs;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  if (vc->hc->type == VRNA_HC_WINDOW)
    return E_stack_window(vc, i, j);

  cp                = vc->cutpoint;
  P                 = vc->params;
  md                = &(P->model_details);
  sn                = vc->strand_number;
  rtype             = &(md->rtype[0]);
  indx              = vc->jindx;
  hard_constraints  = vc->hc->matrix;

  if (vc->hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = vc->hc->f;
    hc_dat_local.hc_dat = vc->hc->data;
  } else {
    evaluate = &hc_default;
  }

  e         = INF;
  p         = i + 1;
  q         = j - 1;
  ij        = indx[j] + i;
  pq        = indx[q] + p;
  eval_loop = (hard_constraints[pq] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) &&
              (hard_constraints[ij] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP);

  if ((j - i - 1) < 2)
    return e;

  if (eval_loop && evaluate(i, j, p, q, VRNA_DECOMP_PAIR_IL, &hc_dat_local)) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        S       = vc->sequence_encoding;
        ptype   = vc->ptype;
        type    = (unsigned char)ptype[ij];
        type_2  = rtype[(unsigned char)ptype[pq]];
        sc      = vc->sc;

        if (type == 0)
          type = 7;

        if (type_2 == 0)
          type_2 = 7;

        if ((sn[p] == sn[i]) && (sn[j] == sn[q])) {
          /* regular stack */
          e = P->stack[type][type_2];
        } else {
          /* stack like cofold structure */
          short si, sj;
          si  = (sn[i + 1] == sn[i]) ? S[i + 1] : -1;
          sj  = (sn[j] == sn[j - 1]) ? S[j - 1] : -1;
          e   = E_IntLoop_Co(rtype[type], rtype[type_2],
                             i, j, p, q,
                             cp,
                             si, sj,
                             S[p - 1], S[q + 1],
                             md->dangles,
                             P);
        }

        /* add soft constraints */
        if (sc) {
          if (sc->energy_bp)
            e += sc->energy_bp[ij];

          if (sc->energy_stack) {
            e += sc->energy_stack[i] +
                 sc->energy_stack[p] +
                 sc->energy_stack[q] +
                 sc->energy_stack[j];
          }

          if (sc->f)
            e += sc->f(i, j, p, q, VRNA_DECOMP_PAIR_IL, sc->data);
        }

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        n_seq = vc->n_seq;
        SS    = vc->S;
        a2s   = vc->a2s;
        scs   = vc->scs;
        e     = 0;
        for (s = 0; s < n_seq; s++) {
          type    = get_pair_type(SS[s][i], SS[s][j], md);
          type_2  = get_pair_type(SS[s][q], SS[s][p], md);  /* q,p not p,q! */
          e       += P->stack[type][type_2];
        }

        if (scs) {
          for (s = 0; s < n_seq; s++) {
            if (scs[s]) {
              if (scs[s]->energy_bp)
                e += scs[s]->energy_bp[ij];

              if (scs[s]->energy_stack) {
                if (SS[s][i] && SS[s][j] && SS[s][p] && SS[s][q]) {
                  /* don't allow gaps in stack */
                  e += scs[s]->energy_stack[a2s[s][i]] +
                       scs[s]->energy_stack[a2s[s][p]] +
                       scs[s]->energy_stack[a2s[s][q]] +
                       scs[s]->energy_stack[a2s[s][j]];
                }
              }

              if (scs[s]->f) {
                e +=
                  scs[s]->f(a2s[s][i],
                            a2s[s][j],
                            a2s[s][p],
                            a2s[s][q],
                            VRNA_DECOMP_PAIR_IL,
                            scs[s]->data);
              }
            }
          }
        }

        break;

      default:
        break;
    }
  }

  return e;
}


PRIVATE int
E_stack_window(vrna_fold_compound_t *vc,
               int                  i,
               int                  j)
{
  char                      **ptype;
  unsigned char             **hard_constraints, eval_loop;
  unsigned int              **a2s;
  short                     **SS;
  int                       e, p, q, *rtype, type, type_2, s, n_seq;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_sc_t                 *sc, **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  P                 = vc->params;
  md                = &(P->model_details);
  rtype             = &(md->rtype[0]);
  hard_constraints  = vc->hc->matrix_local;

  if (vc->hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = vc->hc->f;
    hc_dat_local.hc_dat = vc->hc->data;
  } else {
    evaluate = &hc_default;
  }

  e         = INF;
  p         = i + 1;
  q         = j - 1;
  eval_loop = (hard_constraints[p][q - p] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC) &&
              (hard_constraints[i][j - i] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP);

  if ((j - i - 1) < 2)
    return e;

  if (eval_loop && evaluate(i, j, p, q, VRNA_DECOMP_PAIR_IL, &hc_dat_local)) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        ptype   = vc->ptype_local;
        type    = ptype[i][j - i];
        type_2  = rtype[ptype[p][q - p]];
        sc      = vc->sc;

        if (type == 0)
          type = 7;

        if (type_2 == 0)
          type_2 = 7;

        /* regular stack */
        e = P->stack[type][type_2];

        /* add soft constraints */
        if (sc) {
          if (sc->energy_stack) {
            e += sc->energy_stack[i] +
                 sc->energy_stack[p] +
                 sc->energy_stack[q] +
                 sc->energy_stack[j];
          }

          if (sc->f)
            e += sc->f(i, j, p, q, VRNA_DECOMP_PAIR_IL, sc->data);
        }

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        n_seq = vc->n_seq;
        SS    = vc->S;
        a2s   = vc->a2s;
        scs   = vc->scs;
        e     = 0;
        for (s = 0; s < n_seq; s++) {
          type    = get_pair_type(SS[s][i], SS[s][j], md);
          type_2  = get_pair_type(SS[s][q], SS[s][p], md);  /* q,p not p,q! */
          e       += P->stack[type][type_2];
        }

        /* add soft constraints */
        if (scs) {
          for (s = 0; s < n_seq; s++) {
            if (scs[s]) {
              if (scs[s]->energy_stack) {
                if (SS[s][i] && SS[s][j] && SS[s][p] && SS[s][q]) {
                  /* don't allow gaps in stack */
                  e += scs[s]->energy_stack[a2s[s][i]] +
                       scs[s]->energy_stack[a2s[s][p]] +
                       scs[s]->energy_stack[a2s[s][q]] +
                       scs[s]->energy_stack[a2s[s][j]];
                }
              }

              if (scs[s]->f) {
                e +=
                  scs[s]->f(a2s[s][i],
                            a2s[s][j],
                            a2s[s][p],
                            a2s[s][q],
                            VRNA_DECOMP_PAIR_IL,
                            scs[s]->data);
              }
            }
          }
        }

        break;
    }
  }

  return e;
}


PUBLIC int
vrna_BT_stack(vrna_fold_compound_t  *vc,
              int                   *i,
              int                   *j,
              int                   *en,
              vrna_bp_stack_t       *bp_stack,
              int                   *stack_count)
{
  if (vc) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          return BT_stack_window(vc, i, j, en, bp_stack, stack_count);
        else
          return BT_stack(vc, i, j, en, bp_stack, stack_count);

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          return BT_stack_window_comparative(vc, i, j, en, bp_stack, stack_count);
        else
          return BT_stack_comparative(vc, i, j, en, bp_stack, stack_count);

        break;
    }
  }

  return 0;
}


PRIVATE int
BT_stack(vrna_fold_compound_t *vc,
         int                  *i,
         int                  *j,
         int                  *en,
         vrna_bp_stack_t      *bp_stack,
         int                  *stack_count)
{
  unsigned char             type, type_2;
  char                      *ptype;
  unsigned char             eval_loop;
  unsigned int              *sn;
  int                       ij, p, q, *idx, *my_c, *rtype, cp;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  cp    = vc->cutpoint;
  idx   = vc->jindx;
  P     = vc->params;
  md    = &(P->model_details);
  hc    = vc->hc;
  sc    = vc->sc;
  sn    = vc->strand_number;
  my_c  = vc->matrices->c;
  ij    = idx[*j] + *i;
  ptype = vc->ptype;
  type  = (unsigned char)ptype[ij];
  rtype = &(md->rtype[0]);
  p     = *i + 1;
  q     = *j - 1;

  if (vc->hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = vc->hc->f;
    hc_dat_local.hc_dat = vc->hc->data;
  } else {
    evaluate = &hc_default;
  }

  if (my_c[ij] == *en) {
    /*  always true, if (i.j) closes canonical structure,
     * thus (i+1.j-1) must be a pair
     */
    eval_loop = (hc->matrix[ij] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP)
                && (hc->matrix[idx[q] + p] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC);

    if (eval_loop && evaluate(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, &hc_dat_local)) {
      type_2  = ptype[idx[q] + p];
      type_2  = rtype[type_2];

      if (type == 0)
        type = 7;

      if (type_2 == 0)
        type_2 = 7;

      if ((sn[p] == sn[*i]) && (sn[*j] == sn[q])) {
        /* regular stack */
        *en -= P->stack[type][type_2];
      } else {
        /* stack like cofold structure */
        short si, sj, *S;
        S   = vc->sequence_encoding;
        si  = (sn[p] == sn[*i]) ? S[p] : -1;
        sj  = (sn[*j] == sn[q]) ? S[q] : -1;
        *en -= E_IntLoop_Co(rtype[type], rtype[type_2],
                            *i, *j, p, q,
                            cp,
                            si, sj,
                            S[p - 1], S[q + 1],
                            md->dangles,
                            P);
      }

      if (sc) {
        if (sc->energy_bp)
          *en -= sc->energy_bp[ij];

        if (sc->energy_stack) {
          *en -= sc->energy_stack[*i] +
                 sc->energy_stack[p] +
                 sc->energy_stack[q] +
                 sc->energy_stack[*j];
        }

        if (sc->f)
          *en -= sc->f(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, sc->data);
      }

      bp_stack[++(*stack_count)].i  = p;
      bp_stack[(*stack_count)].j    = q;
      (*i)++;
      (*j)--;
      return 1;
    }
  }

  return 0;
}


PRIVATE int
BT_stack_comparative(vrna_fold_compound_t *vc,
                     int                  *i,
                     int                  *j,
                     int                  *en,
                     vrna_bp_stack_t      *bp_stack,
                     int                  *stack_count)
{
  short                     **S;
  int                       type, type_2;
  unsigned char             eval_loop;
  int                       p, q, *c, n_seq, ss, ij, *idx;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n_seq = vc->n_seq;
  S     = vc->S;
  P     = vc->params;
  md    = &(P->model_details);
  hc    = vc->hc;
  scs   = vc->scs;
  c     = vc->matrices->c;
  idx   = vc->jindx;
  ij    = idx[*j] + *i;
  p     = *i + 1;
  q     = *j - 1;

  if (vc->hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = vc->hc->f;
    hc_dat_local.hc_dat = vc->hc->data;
  } else {
    evaluate = &hc_default;
  }

  if (c[ij] == *en) {
    /*  always true, if (i.j) closes canonical structure,
     * thus (i+1.j-1) must be a pair
     */
    eval_loop = (hc->matrix[ij] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP)
                && (hc->matrix[idx[q] + p] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC);

    if (eval_loop && evaluate(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, &hc_dat_local)) {
      for (ss = 0; ss < n_seq; ss++) {
        type    = get_pair_type(S[ss][*i], S[ss][*j], md);
        type_2  = get_pair_type(S[ss][q], S[ss][p], md);
        *en     -= P->stack[type][type_2];
      }

      if (scs) {
        for (ss = 0; ss < n_seq; ss++)
          if (scs[ss]) {
            if (scs[ss]->energy_stack) {
              *en -= scs[ss]->energy_stack[*i] +
                     scs[ss]->energy_stack[p] +
                     scs[ss]->energy_stack[q] +
                     scs[ss]->energy_stack[*j];
            }

            if (scs[ss]->f)
              *en -= scs[ss]->f(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, scs[ss]->data);
          }
      }

      *en += vc->pscore[ij];

      bp_stack[++(*stack_count)].i  = p;
      bp_stack[(*stack_count)].j    = q;
      (*i)++;
      (*j)--;
      return 1;
    }
  }

  return 0;
}


PRIVATE int
BT_stack_window(vrna_fold_compound_t  *vc,
                int                   *i,
                int                   *j,
                int                   *en,
                vrna_bp_stack_t       *bp_stack,
                int                   *stack_count)
{
  unsigned char             type, type_2;
  char                      **ptype;
  unsigned char             eval_loop;
  int                       p, q, **c, *rtype;
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
  c     = vc->matrices->c_local;
  ptype = vc->ptype_local;
  type  = (unsigned char)ptype[*i][*j - *i];
  rtype = &(md->rtype[0]);
  p     = *i + 1;
  q     = *j - 1;

  if (vc->hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = vc->hc->f;
    hc_dat_local.hc_dat = vc->hc->data;
  } else {
    evaluate = &hc_default;
  }

  if (c[*i][*j - *i] == *en) {
    /*  always true, if (i.j) closes canonical structure,
     * thus (i+1.j-1) must be a pair
     */
    eval_loop = (hc->matrix_local[*i][*j - *i] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP)
                && (hc->matrix_local[p][q - p] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC);

    if (eval_loop && evaluate(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, &hc_dat_local)) {
      type_2  = ptype[p][q - p];
      type_2  = rtype[type_2];

      if (type == 0)
        type = 7;

      if (type_2 == 0)
        type_2 = 7;

      *en -= P->stack[type][type_2];

      if (sc) {
        if (sc->energy_stack) {
          *en -= sc->energy_stack[*i] +
                 sc->energy_stack[p] +
                 sc->energy_stack[q] +
                 sc->energy_stack[*j];
        }

        if (sc->f)
          *en -= sc->f(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, sc->data);
      }

      bp_stack[++(*stack_count)].i  = p;
      bp_stack[(*stack_count)].j    = q;
      (*i)++;
      (*j)--;
      return 1;
    }
  }

  return 0;
}


PRIVATE int
BT_stack_window_comparative(vrna_fold_compound_t  *vc,
                            int                   *i,
                            int                   *j,
                            int                   *en,
                            vrna_bp_stack_t       *bp_stack,
                            int                   *stack_count)
{
  int                       type, type_2;
  unsigned char             eval_loop;
  short                     **S;
  int                       p, q, **c, n_seq, ss;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n_seq = vc->n_seq;
  S     = vc->S;
  P     = vc->params;
  md    = &(P->model_details);
  hc    = vc->hc;
  scs   = vc->scs;
  c     = vc->matrices->c_local;
  p     = *i + 1;
  q     = *j - 1;

  if (vc->hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = vc->hc->f;
    hc_dat_local.hc_dat = vc->hc->data;
  } else {
    evaluate = &hc_default;
  }

  if (c[*i][*j - *i] == *en) {
    /*  always true, if (i.j) closes canonical structure,
     * thus (i+1.j-1) must be a pair
     */
    eval_loop = (hc->matrix_local[*i][*j - *i] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP)
                && (hc->matrix_local[p][q - p] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC);

    if (eval_loop && evaluate(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, &hc_dat_local)) {
      for (ss = 0; ss < n_seq; ss++) {
        type    = get_pair_type(S[ss][*i], S[ss][*j], md);
        type_2  = get_pair_type(S[ss][q], S[ss][p], md);
        *en     -= P->stack[type][type_2];
      }

      if (scs) {
        for (ss = 0; ss < n_seq; ss++)
          if (scs[ss]) {
            if (scs[ss]->energy_stack) {
              *en -= scs[ss]->energy_stack[*i] +
                     scs[ss]->energy_stack[p] +
                     scs[ss]->energy_stack[q] +
                     scs[ss]->energy_stack[*j];
            }

            if (scs[ss]->f)
              *en -= scs[ss]->f(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, scs[ss]->data);
          }
      }

      *en += vc->pscore_local[*i][*j - *i];

      bp_stack[++(*stack_count)].i  = p;
      bp_stack[(*stack_count)].j    = q;
      (*i)++;
      (*j)--;
      return 1;
    }
  }

  return 0;
}


PUBLIC int
vrna_BT_int_loop(vrna_fold_compound_t *vc,
                 int                  *i,
                 int                  *j,
                 int                  en,
                 vrna_bp_stack_t      *bp_stack,
                 int                  *stack_count)
{
  if (vc) {
    switch (vc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          return BT_int_loop_window(vc, i, j, en, bp_stack, stack_count);
        else
          return BT_int_loop(vc, i, j, en, bp_stack, stack_count);

        break;

      case VRNA_FC_TYPE_COMPARATIVE:
        if (vc->hc->type == VRNA_HC_WINDOW)
          return BT_int_loop_window_comparative(vc, i, j, en, bp_stack, stack_count);
        else
          return BT_int_loop_comparative(vc, i, j, en, bp_stack, stack_count);

        break;
    }
  }

  return 0;
}


PRIVATE int
BT_int_loop(vrna_fold_compound_t  *vc,
            int                   *i,
            int                   *j,
            int                   en,
            vrna_bp_stack_t       *bp_stack,
            int                   *stack_count)
{
  unsigned char             type, type_2;
  char                      *ptype;
  unsigned char             eval_loop;
  short                     *S1;
  unsigned int              *sn;
  int                       ij, p, q, minq, turn, *idx, noGUclosure, no_close,
                            energy, new, *my_c, *rtype;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_ud_t                 *domains_up;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  idx         = vc->jindx;
  P           = vc->params;
  md          = &(P->model_details);
  hc          = vc->hc;
  sc          = vc->sc;
  sn          = vc->strand_number;
  my_c        = vc->matrices->c;
  turn        = md->min_loop_size;
  ij          = idx[*j] + *i;
  ptype       = vc->ptype;
  type        = (unsigned char)ptype[ij];
  rtype       = &(md->rtype[0]);
  S1          = vc->sequence_encoding;
  noGUclosure = md->noGUclosure;
  no_close    = (((type == 3) || (type == 4)) && noGUclosure);
  domains_up  = vc->domains_up;

  if (vc->hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = vc->hc->f;
    hc_dat_local.hc_dat = vc->hc->data;
  } else {
    evaluate = &hc_default;
  }

  if (hc->matrix[ij] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
    if (type == 0)
      type = 7;

    if (domains_up && domains_up->energy_cb) {
      for (p = *i + 1; p <= MIN2(*j - 2 - turn, *i + MAXLOOP + 1); p++) {
        minq = *j - *i + p - MAXLOOP - 2;
        if (minq < p + 1 + turn)
          minq = p + 1 + turn;

        if (hc->up_int[*i + 1] < (p - *i - 1))
          break;

        for (q = *j - 1; q >= minq; q--) {
          if (hc->up_int[q + 1] < (*j - q - 1))
            break;

          type_2    = (unsigned char)ptype[idx[q] + p];
          eval_loop = hc->matrix[idx[q] + p] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC;

          if (!(eval_loop && evaluate(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, &hc_dat_local)))
            continue;

          type_2 = rtype[type_2];

          if (type_2 == 0)
            type_2 = 7;

          if (noGUclosure)
            if (no_close || (type_2 == 3) || (type_2 == 4))
              if ((p > *i + 1) || (q < *j - 1))
                continue;

          /* continue unless stack */

          energy  = eval_interior_loop(vc, *i, *j, p, q);
          new     = energy + my_c[idx[q] + p];

          if (new == en) {
            bp_stack[++(*stack_count)].i  = p;
            bp_stack[(*stack_count)].j    = q;
            if (sc) {
              if (sc->bt) {
                vrna_basepair_t *ptr, *aux_bps;
                aux_bps = sc->bt(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, sc->data);
                for (ptr = aux_bps; ptr && ptr->i != 0; ptr++) {
                  bp_stack[++(*stack_count)].i  = ptr->i;
                  bp_stack[(*stack_count)].j    = ptr->j;
                }
                free(aux_bps);
              }
            }

            *i = p, *j = q;
            return 1; /* success */
          }
        }
      }
    } else {
      for (p = *i + 1; p <= MIN2(*j - 2 - turn, *i + MAXLOOP + 1); p++) {
        minq = *j - *i + p - MAXLOOP - 2;
        if (minq < p + 1 + turn)
          minq = p + 1 + turn;

        if (hc->up_int[*i + 1] < (p - *i - 1))
          break;

        for (q = *j - 1; q >= minq; q--) {
          if (hc->up_int[q + 1] < (*j - q - 1))
            break;

          type_2 = (unsigned char)ptype[idx[q] + p];

          eval_loop = hc->matrix[idx[q] + p] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC;

          if (!(eval_loop && evaluate(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, &hc_dat_local)))
            continue;

          type_2 = rtype[type_2];

          if (type_2 == 0)
            type_2 = 7;

          if (noGUclosure)
            if (no_close || (type_2 == 3) || (type_2 == 4))
              if ((p > *i + 1) || (q < *j - 1))
                continue;

          /* continue unless stack */

          energy = ubf_eval_int_loop(*i, *j, p, q,
                                     (*i) + 1, (*j) - 1, p - 1, q + 1,
                                     S1[*i + 1], S1[*j - 1], S1[p - 1], S1[q + 1],
                                     type, type_2,
                                     rtype,
                                     ij,
                                     -1,
                                     P,
                                     sc);
          new = energy + my_c[idx[q] + p];

          if (new == en) {
            bp_stack[++(*stack_count)].i  = p;
            bp_stack[(*stack_count)].j    = q;
            if (sc) {
              if (sc->bt) {
                vrna_basepair_t *ptr, *aux_bps;
                aux_bps = sc->bt(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, sc->data);
                for (ptr = aux_bps; ptr && ptr->i != 0; ptr++) {
                  bp_stack[++(*stack_count)].i  = ptr->i;
                  bp_stack[(*stack_count)].j    = ptr->j;
                }
                free(aux_bps);
              }
            }

            *i = p, *j = q;
            return 1; /* success */
          }
        }
      }
    }
  }

  /* is it a g-quadruplex? */
  if (md->gquad) {
    /*
     * The case that is handled here actually resembles something like
     * an interior loop where the enclosing base pair is of regular
     * kind and the enclosed pair is not a canonical one but a g-quadruplex
     * that should then be decomposed further...
     */
    if (sn[*j] == sn[*i]) {
      if (vrna_BT_gquad_int(vc, *i, *j, en, bp_stack, stack_count)) {
        *i = *j = -1; /* tell the calling block to continue backtracking with next block */
        return 1;
      }
    }
  }

  return 0; /* unsuccessful */
}


PRIVATE int
BT_int_loop_comparative(vrna_fold_compound_t  *vc,
                        int                   *i,
                        int                   *j,
                        int                   en,
                        vrna_bp_stack_t       *bp_stack,
                        int                   *stack_count)
{
  unsigned char             eval_loop;
  unsigned int              **a2s;
  short                     **S, **S5, **S3, *S_cons;
  int                       ij, p, q, minq, turn, *idx,
                            energy, new, *my_c, *ggg, n_seq, ss, *type, type_2;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 **scs;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n_seq   = vc->n_seq;
  S_cons  = vc->S_cons;
  S       = vc->S;
  S5      = vc->S5;
  S3      = vc->S3;
  a2s     = vc->a2s;
  idx     = vc->jindx;
  P       = vc->params;
  md      = &(P->model_details);
  hc      = vc->hc;
  scs     = vc->scs;
  my_c    = vc->matrices->c;
  ggg     = vc->matrices->ggg;
  turn    = md->min_loop_size;
  ij      = idx[*j] + *i;

  if (vc->hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = vc->hc->f;
    hc_dat_local.hc_dat = vc->hc->data;
  } else {
    evaluate = &hc_default;
  }

  if (hc->matrix[ij] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
    type = (int *)vrna_alloc(n_seq * sizeof(int));
    for (ss = 0; ss < n_seq; ss++)
      type[ss] = get_pair_type(S[ss][*i], S[ss][*j], md);

    for (p = *i + 1; p <= MIN2(*j - 2 - turn, *i + MAXLOOP + 1); p++) {
      minq = *j - *i + p - MAXLOOP - 2;
      if (minq < p + 1 + turn)
        minq = p + 1 + turn;

      if (hc->up_int[*i + 1] < (p - *i - 1))
        break;

      for (q = *j - 1; q >= minq; q--) {
        if (hc->up_int[q + 1] < (*j - q - 1))
          break;

        eval_loop = hc->matrix[idx[q] + p] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC;

        if (!(eval_loop && evaluate(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, &hc_dat_local)))
          continue;

        for (ss = energy = 0; ss < n_seq; ss++) {
          int u1  = a2s[ss][p - 1] - a2s[ss][*i];
          int u2  = a2s[ss][*j - 1] - a2s[ss][q];
          type_2  = get_pair_type(S[ss][q], S[ss][p], md); /* q,p not p,q */
          energy  += E_IntLoop(u1,
                               u2,
                               type[ss],
                               type_2,
                               S3[ss][*i],
                               S5[ss][*j],
                               S5[ss][p],
                               S3[ss][q],
                               P);
        }

        if (scs) {
          for (ss = 0; ss < n_seq; ss++) {
            if (scs[ss]) {
              int u1  = a2s[ss][p - 1] - a2s[ss][*i];
              int u2  = a2s[ss][*j - 1] - a2s[ss][q];
              /*
               *            int u1 = p - i - 1;
               *            int u2 = j - q - 1;
               */
              if (u1 + u2 == 0) {
                if (scs[ss]->energy_stack) {
                  if (S[ss][*i] && S[ss][*j] && S[ss][p] && S[ss][q]) {
                    /* don't allow gaps in stack */
                    energy += scs[ss]->energy_stack[a2s[ss][*i]] +
                              scs[ss]->energy_stack[a2s[ss][p]] +
                              scs[ss]->energy_stack[a2s[ss][q]] +
                              scs[ss]->energy_stack[a2s[ss][*j]];
                  }
                }
              }

              if (scs[ss]->energy_bp)
                energy += scs[ss]->energy_bp[ij];

              if (scs[ss]->energy_up)
                energy += scs[ss]->energy_up[a2s[ss][*i] + 1][u1] +
                          scs[ss]->energy_up[a2s[ss][q] + 1][u2];
            }
          }
        }

        new = energy + my_c[idx[q] + p];

        if (new == en) {
          bp_stack[++(*stack_count)].i  = p;
          bp_stack[(*stack_count)].j    = q;
          if (scs && scs[0]) {
            if (scs[0]->bt) {
              vrna_basepair_t *ptr, *aux_bps;
              aux_bps = scs[0]->bt(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, scs[0]->data);
              for (ptr = aux_bps; ptr && ptr->i != 0; ptr++) {
                bp_stack[++(*stack_count)].i  = ptr->i;
                bp_stack[(*stack_count)].j    = ptr->j;
              }
              free(aux_bps);
            }
          }

          free(type);
          *i = p, *j = q;
          return 1; /* success */
        }
      }
    }

    free(type);
  }

  /* is it a g-quadruplex? */
  if (md->gquad) {
    /*
     * The case that is handled here actually resembles something like
     * an interior loop where the enclosing base pair is of regular
     * kind and the enclosed pair is not a canonical one but a g-quadruplex
     * that should then be decomposed further...
     */
    type = (int *)vrna_alloc(n_seq * sizeof(int));
    for (ss = 0; ss < n_seq; ss++)
      type[ss] = get_pair_type(S[ss][*i], S[ss][*j], md);

    if (backtrack_GQuad_IntLoop_comparative(en, *i, *j, type, S_cons, S5, S3, ggg, idx, &p, &q,
                                            n_seq,
                                            P)) {
      if (vrna_BT_gquad_mfe(vc, p, q, bp_stack, stack_count)) {
        *i = *j = -1; /* tell the calling block to continue backtracking with next block */
        return 1;
      }
    }

    free(type);
  }

  return 0; /* unsuccessful */
}


PRIVATE int
BT_int_loop_window(vrna_fold_compound_t *vc,
                   int                  *i,
                   int                  *j,
                   int                  en,
                   vrna_bp_stack_t      *bp_stack,
                   int                  *stack_count)
{
  int                       type, type_2;
  char                      **ptype;
  unsigned char             eval_loop;
  short                     *S, *S1;
  unsigned int              *sn;
  int                       p, q, minq, turn, maxdist, noGUclosure, no_close,
                            energy, new, **c, **ggg, *rtype, u1, u2;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  S1          = vc->sequence_encoding;
  S           = vc->sequence_encoding2;
  sn          = vc->strand_number;
  maxdist     = vc->window_size;
  ptype       = vc->ptype_local;
  P           = vc->params;
  md          = &(P->model_details);
  hc          = vc->hc;
  sc          = vc->sc;
  c           = vc->matrices->c_local;
  ggg         = vc->matrices->ggg_local;
  turn        = md->min_loop_size;
  noGUclosure = md->noGUclosure;
  rtype       = &(md->rtype[0]);
  type        = ptype[*i][*j - *i];
  no_close    = (((type == 3) || (type == 4)) && noGUclosure);

  if (hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default;
  }

  if (hc->matrix_local[*i][*j - *i] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
    if (type == 0)
      type = 7;

    for (p = *i + 1; p <= MIN2(*j - 2 - turn, *i + MAXLOOP + 1); p++) {
      u1 = p - (*i) - 1;

      minq = *j - *i + p - MAXLOOP - 2;
      if (minq < p + 1 + turn)
        minq = p + 1 + turn;

      if (hc->up_int[*i + 1] < u1)
        break;

      for (q = *j - 1; q >= minq; q--) {
        u2 = *j - q - 1;

        if (hc->up_int[q + 1] < u2)
          break;

        eval_loop = hc->matrix_local[p][q - p] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC;

        if (!(eval_loop && evaluate(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, &hc_dat_local)))
          continue;

        type_2  = ptype[p][q - p];
        type_2  = rtype[type_2];

        if (type_2 == 0)
          type_2 = 7;

        if (noGUclosure)
          if (no_close || (type_2 == 3) || (type_2 == 4))
            if ((p > *i + 1) || (q < *j - 1))
              continue;

        /* continue unless stack */

        energy = E_IntLoop(u1,
                           u2,
                           type,
                           type_2,
                           S1[*i + 1],
                           S1[*j - 1],
                           S1[p - 1],
                           S1[q + 1],
                           P);

        if (sc) {
          if (sc->energy_up)
            energy += sc->energy_up[*i + 1][u1] +
                      sc->energy_up[q + 1][u2];

          if (sc->energy_bp_local)
            energy += sc->energy_bp_local[*i][*j - *i];

          if (sc->energy_stack) {
            if (u1 + u2 == 0) {
              energy += sc->energy_stack[*i] +
                        sc->energy_stack[p] +
                        sc->energy_stack[q] +
                        sc->energy_stack[*j];
            }
          }

          if (sc->f)
            energy += sc->f(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, sc->data);
        }

        new = energy + c[p][q - p];

        if (new == en) {
          bp_stack[++(*stack_count)].i  = p;
          bp_stack[(*stack_count)].j    = q;
          if (sc) {
            if (sc->bt) {
              vrna_basepair_t *ptr, *aux_bps;
              aux_bps = sc->bt(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, sc->data);
              for (ptr = aux_bps; ptr && ptr->i != 0; ptr++) {
                bp_stack[++(*stack_count)].i  = ptr->i;
                bp_stack[(*stack_count)].j    = ptr->j;
              }
              free(aux_bps);
            }
          }

          *i = p, *j = q;
          return 1; /* success */
        }
      }
    }
  }

  /* is it a g-quadruplex? */
  if (md->gquad) {
    /*
     * The case that is handled here actually resembles something like
     * an interior loop where the enclosing base pair is of regular
     * kind and the enclosed pair is not a canonical one but a g-quadruplex
     * that should then be decomposed further...
     */
    if (sn[*j] == sn[*i]) {
      if (backtrack_GQuad_IntLoop_L(en, *i, *j, type, S, ggg, maxdist, &p, &q, P)) {
        if (vrna_BT_gquad_mfe(vc, p, q, bp_stack, stack_count)) {
          *i = *j = -1; /* tell the calling block to continue backtracking with next block */
          return 1;
        }
      }
    }
  }

  return 0; /* unsuccessful */
}


PRIVATE int
BT_int_loop_window_comparative(vrna_fold_compound_t *vc,
                               int                  *i,
                               int                  *j,
                               int                  en,
                               vrna_bp_stack_t      *bp_stack,
                               int                  *stack_count)
{
  unsigned char             eval_loop;
  unsigned int              **a2s;
  short                     **S, **SS, **S5, **S3, *S_cons;
  int                       *type, type_2, *rtype;
  int                       p, q, minq, turn, energy, new, **c, **ggg, n_seq, ss;
  vrna_param_t              *P;
  vrna_md_t                 *md;
  vrna_hc_t                 *hc;
  vrna_sc_t                 **scs, *sc;
  vrna_callback_hc_evaluate *evaluate;
  struct default_data       hc_dat_local;

  n_seq   = vc->n_seq;
  S_cons  = vc->S_cons;
  S       = vc->S;
  S5      = vc->S5;       /* S5[s][i] holds next base 5' of i in sequence s */
  S3      = vc->S3;       /* Sl[s][i] holds next base 3' of i in sequence s */
  a2s     = vc->a2s;
  P       = vc->params;
  md      = &(P->model_details);
  hc      = vc->hc;
  scs     = vc->scs;
  c       = vc->matrices->c_local;
  ggg     = vc->matrices->ggg_local;
  turn    = md->min_loop_size;
  rtype   = &(md->rtype[0]);

  if (hc->f) {
    evaluate            = &hc_default_user;
    hc_dat_local.hc_f   = hc->f;
    hc_dat_local.hc_dat = hc->data;
  } else {
    evaluate = &hc_default;
  }

  if (hc->matrix_local[*i][*j - *i] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
    type = (int *)vrna_alloc(n_seq * sizeof(int));
    for (ss = 0; ss < n_seq; ss++)
      type[ss] = get_pair_type(S[ss][*i], S[ss][*j], md);

    for (p = *i + 1; p <= MIN2(*j - 2 - turn, *i + MAXLOOP + 1); p++) {
      minq = *j - *i + p - MAXLOOP - 2;
      if (minq < p + 1 + turn)
        minq = p + 1 + turn;

      if (hc->up_int[*i + 1] < (p - *i - 1))
        break;

      for (q = *j - 1; q >= minq; q--) {
        if (hc->up_int[q + 1] < (*j - q - 1))
          break;

        eval_loop = hc->matrix_local[p][q - p] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP_ENC;

        if (!(eval_loop && evaluate(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, &hc_dat_local)))
          continue;

        for (ss = energy = 0; ss < n_seq; ss++) {
          type_2  = get_pair_type(S[ss][q], S[ss][p], md); /* q,p not p,q! */
          sc      = (scs && scs[ss]) ? scs[ss] : NULL;
          energy  += ubf_eval_int_loop_comparative(*i, *j, p, q,
                                                   type[ss],
                                                   type_2,
                                                   rtype,
                                                   0, -1,
                                                   P,
                                                   S[ss],
                                                   S5[ss],
                                                   S3[ss],
                                                   a2s[ss],
                                                   sc);
        }

        new = energy + c[p][q - p];

        if (new == en) {
          bp_stack[++(*stack_count)].i  = p;
          bp_stack[(*stack_count)].j    = q;
          if (scs && scs[0]) {
            if (scs[0]->bt) {
              vrna_basepair_t *ptr, *aux_bps;
              aux_bps = scs[0]->bt(*i, *j, p, q, VRNA_DECOMP_PAIR_IL, scs[0]->data);
              for (ptr = aux_bps; ptr && ptr->i != 0; ptr++) {
                bp_stack[++(*stack_count)].i  = ptr->i;
                bp_stack[(*stack_count)].j    = ptr->j;
              }
              free(aux_bps);
            }
          }

          free(type);
          *i = p, *j = q;
          return 1; /* success */
        }
      }
    }

    free(type);
  }

  /* is it a g-quadruplex? */
  if (md->gquad) {
    /*
     * The case that is handled here actually resembles something like
     * an interior loop where the enclosing base pair is of regular
     * kind and the enclosed pair is not a canonical one but a g-quadruplex
     * that should then be decomposed further...
     */
    type = (int *)vrna_alloc(n_seq * sizeof(int));
    for (ss = 0; ss < n_seq; ss++)
      type[ss] = get_pair_type(S[ss][*i], S[ss][*j], md);

    if (backtrack_GQuad_IntLoop_L_comparative(en, *i, *j, type, S_cons, S5, S3, ggg, &p, &q, n_seq,
                                              P)) {
      if (vrna_BT_gquad_mfe(vc, p, q, bp_stack, stack_count)) {
        *i = *j = -1; /* tell the calling block to continue backtracking with next block */
        return 1;
      }
    }

    free(type);
  }

  return 0; /* unsuccessful */
}


PRIVATE unsigned char
hc_default(int            i,
           int            j,
           int            k,
           int            l,
           unsigned char  d,
           void           *data)
{
  return (unsigned char)1;
}


PRIVATE unsigned char
hc_default_user(int           i,
                int           j,
                int           k,
                int           l,
                unsigned char d,
                void          *data)
{
  struct default_data *dat = (struct default_data *)data;

  return dat->hc_f(i, j, k, l, d, dat->hc_dat);
}
