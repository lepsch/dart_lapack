      void slarre(final int RANGE, final int N, final int VL, final int VU, final int IL, final int IU, final int D, final int E, final int E2, final int RTOL1, final int RTOL2, final int SPLTOL, final int NSPLIT, final int ISPLIT, final int M, final int W, final int WERR, final int WGAP, final int IBLOCK, final int INDEXW, final int GERS, final int PIVMIN, final Array<double> _WORK_, final Array<int> IWORK_, final Box<int> INFO,) {
  final _WORK = _WORK_.dim();
  final IWORK = IWORK_.dim();

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             RANGE;
      int                IL, INFO, IU, M, N, NSPLIT;
      double               PIVMIN, RTOL1, RTOL2, SPLTOL, VL, VU;
      int                IBLOCK( * ), ISPLIT( * ), IWORK( * ), INDEXW( * );
      double               D( * ), E( * ), E2( * ), GERS( * ), W( * ),WERR( * ), WGAP( * ), WORK( * );
      // ..

      double               FAC, FOUR, FOURTH, FUDGE, HALF, HNDRD, MAXGROWTH, ONE, PERT, TWO, ZERO       PARAMETER          ( ZERO = 0.0, ONE = 1.0, TWO = 2.0, FOUR=4.0, HNDRD = 100.0, PERT = 4.0, HALF = ONE/TWO, FOURTH = ONE/FOUR, FAC= HALF, MAXGROWTH = 64.0, FUDGE = 2.0 );
      int                MAXTRY, ALLRNG, INDRNG, VALRNG;
      const              MAXTRY = 6, ALLRNG = 1, INDRNG = 2, VALRNG = 3 ;
      bool               FORCEB, NOREP, USEDQD;
      int                CNT, CNT1, CNT2, I, IBEGIN, IDUM, IEND, IINFO, IN, INDL, INDU, IRANGE, J, JBLK, MB, MM, WBEGIN, WEND;
      double               AVGAP, BSRTOL, CLWDTH, DMAX, DPIVOT, EABS, EMAX, EOLD, EPS, GL, GU, ISLEFT, ISRGHT, RTL, RTOL, S1, S2, SAFMIN, SGNDEF, SIGMA, SPDIAM, TAU, TMP, TMP1;


      int                ISEED( 4 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL                        SLAMCH;
      // EXTERNAL SLAMCH, lsame

      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SLARNV, SLARRA, SLARRB, SLARRC, SLARRD, SLASQ2, SLARRK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN



      INFO = 0;
      NSPLIT = 0;
      M = 0;

      // Quick return if possible

      if ( N <= 0 ) {
         return;
      }

      // Decode RANGE

      if ( lsame( RANGE, 'A' ) ) {
         IRANGE = ALLRNG;
      } else if ( lsame( RANGE, 'V' ) ) {
         IRANGE = VALRNG;
      } else if ( lsame( RANGE, 'I' ) ) {
         IRANGE = INDRNG;
      }

      // Get machine constants
      SAFMIN = SLAMCH( 'S' );
      EPS = SLAMCH( 'P' );

      // Set parameters
      RTL = HNDRD*EPS;
      // If one were ever to ask for less initial precision in BSRTOL,
      // one should keep in mind that for the subset case, the extremal
      // eigenvalues must be at least as accurate as the current setting
      // (eigenvalues in the middle need not as much accuracy)
      BSRTOL = sqrt(EPS)*(0.5e-3);

      // Treat case of 1x1 matrix for quick return;
      if ( N == 1 ) {
         if ( (IRANGE == ALLRNG) || ((IRANGE == VALRNG) && (D(1) > VL) && (D(1) <= VU)) || ((IRANGE == INDRNG) && (IL == 1) && (IU == 1)) ) {
            M = 1;
            W[1] = D(1);
            // The computation error of the eigenvalue is zero
            WERR[1] = ZERO;
            WGAP[1] = ZERO;
            IBLOCK[1] = 1;
            INDEXW[1] = 1;
            GERS[1] = D( 1 );
            GERS[2] = D( 1 );
         }
         // store the shift for the initial RRR, which is zero in this case
         E[1] = ZERO;
         return;
      }

      // General case: tridiagonal matrix of order > 1

      // Init WERR, WGAP. Compute Gerschgorin intervals and spectral diameter.
      // Compute maximum off-diagonal entry and pivmin.
      GL = D(1);
      GU = D(1);
      EOLD = ZERO;
      EMAX = ZERO;
      E[N] = ZERO;
      for (I = 1; I <= N; I++) { // 5
         WERR[I] = ZERO;
         WGAP[I] = ZERO;
         EABS = ( E(I) ).abs();
         if ( EABS >= EMAX ) {
            EMAX = EABS;
         }
         TMP1 = EABS + EOLD;
         GERS[2*I-1] = D(I) - TMP1;
         GL =  min( GL, GERS( 2*I - 1));
         GERS[2*I] = D(I) + TMP1;
         GU = max( GU, GERS(2*I) );
         EOLD  = EABS;
      } // 5
      // The minimum pivot allowed in the Sturm sequence for T
      PIVMIN = SAFMIN * max( ONE, EMAX**2 );
      // Compute spectral diameter. The Gerschgorin bounds give an
      // estimate that is wrong by at most a factor of sqrt(2)
      SPDIAM = GU - GL;

      // Compute splitting points
      slarra(N, D, E, E2, SPLTOL, SPDIAM, NSPLIT, ISPLIT, IINFO );

      // Can force use of bisection instead of faster DQDS.
      // Option left in the code for future multisection work.
      FORCEB = false;

      // Initialize USEDQD, DQDS should be used for ALLRNG unless someone
      // explicitly wants bisection.
      USEDQD = (( IRANGE == ALLRNG ) && ( !FORCEB));

      if ( (IRANGE == ALLRNG) && ( !FORCEB) ) {
         // Set interval [VL,VU] that contains all eigenvalues
         VL = GL;
         VU = GU;
      } else {
         // We call SLARRD to find crude approximations to the eigenvalues
         // in the desired range. In case IRANGE = INDRNG, we also obtain the
         // interval (VL,VU] that contains all the wanted eigenvalues.
         // An interval [LEFT,RIGHT] has converged if
         // RIGHT-LEFT < RTOL*max(ABS(LEFT),ABS(RIGHT))
         // SLARRD needs a WORK of size 4*N, IWORK of size 3*N
         slarrd(RANGE, 'B', N, VL, VU, IL, IU, GERS, BSRTOL, D, E, E2, PIVMIN, NSPLIT, ISPLIT, MM, W, WERR, VL, VU, IBLOCK, INDEXW, WORK, IWORK, IINFO );
         if ( IINFO != 0 ) {
            INFO = -1;
            return;
         }
         // Make sure that the entries M+1 to N in W, WERR, IBLOCK, INDEXW are 0
         for (I = MM+1; I <= N; I++) { // 14
            W[I] = ZERO;
            WERR[I] = ZERO;
            IBLOCK[I] = 0;
            INDEXW[I] = 0;
         } // 14
      }


// **
      // Loop over unreduced blocks
      IBEGIN = 1;
      WBEGIN = 1;
      for (JBLK = 1; JBLK <= NSPLIT; JBLK++) { // 170
         IEND = ISPLIT( JBLK );
         IN = IEND - IBEGIN + 1;

         // 1 X 1 block
         if ( IN == 1 ) {
            if ( (IRANGE == ALLRNG) || ( (IRANGE == VALRNG) && ( D( IBEGIN ) > VL ) && ( D( IBEGIN ) <= VU ) ) || ( (IRANGE == INDRNG) && (IBLOCK(WBEGIN) == JBLK)) ) {
               M = M + 1;
               W[M] = D( IBEGIN );
               WERR[M] = ZERO;
               // The gap for a single block doesn't matter for the later
               // algorithm and is assigned an arbitrary large value
               WGAP[M] = ZERO;
               IBLOCK[M] = JBLK;
               INDEXW[M] = 1;
               WBEGIN = WBEGIN + 1;
            }
            // E( IEND ) holds the shift for the initial RRR
            E[IEND] = ZERO;
            IBEGIN = IEND + 1;
            GO TO 170;
         }

         // Blocks of size larger than 1x1

         // E( IEND ) will hold the shift for the initial RRR, for now set it =0
         E[IEND] = ZERO;

         // Find local outer bounds GL,GU for the block
         GL = D(IBEGIN);
         GU = D(IBEGIN);
         for (I = IBEGIN; I <= IEND; I++) { // 15
            GL = min( GERS( 2*I-1 ), GL );
            GU = max( GERS( 2*I ), GU );
         } // 15
         SPDIAM = GU - GL;

         if ( !((IRANGE == ALLRNG) && ( !FORCEB)) ) {
            // Count the number of eigenvalues in the current block.
            MB = 0;
            for (I = WBEGIN; I <= MM; I++) { // 20
               if ( IBLOCK(I) == JBLK ) {
                  MB = MB+1;
               } else {
                  GOTO 21;
               }
            } // 20
            } // 21

            if ( MB == 0) {
               // No eigenvalue in the current block lies in the desired range
               // E( IEND ) holds the shift for the initial RRR
               E[IEND] = ZERO;
               IBEGIN = IEND + 1;
               GO TO 170;
            } else {

               // Decide whether dqds or bisection is more efficient
               USEDQD = ( (MB > FAC*IN) && ( !FORCEB) );
               WEND = WBEGIN + MB - 1;
               // Calculate gaps for the current block
               // In later stages, when representations for individual
               // eigenvalues are different, we use SIGMA = E( IEND ).
               SIGMA = ZERO;
               for (I = WBEGIN; I <= WEND - 1; I++) { // 30
                  WGAP[I] = max( ZERO, W(I+1)-WERR(I+1) - (W(I)+WERR(I)) );
               } // 30
               WGAP[WEND] = max( ZERO, VU - SIGMA - (W( WEND )+WERR( WEND )));
               // Find local index of the first and last desired evalue.
               INDL = INDEXW(WBEGIN);
               INDU = INDEXW( WEND );
            }
         }
         if (( (IRANGE == ALLRNG) && ( !FORCEB) ) || USEDQD) {
            // Case of DQDS
            // Find approximations to the extremal eigenvalues of the block
            slarrk(IN, 1, GL, GU, D(IBEGIN), E2(IBEGIN), PIVMIN, RTL, TMP, TMP1, IINFO );
            if ( IINFO != 0 ) {
               INFO = -1;
               return;
            }
            ISLEFT = max(GL, TMP - TMP1 - HNDRD * EPS* (TMP - TMP1).abs());
             slarrk(IN, IN, GL, GU, D(IBEGIN), E2(IBEGIN), PIVMIN, RTL, TMP, TMP1, IINFO );
            if ( IINFO != 0 ) {
               INFO = -1;
               return;
            }
            ISRGHT = min(GU, TMP + TMP1 + HNDRD * EPS * (TMP + TMP1).abs());
            // Improve the estimate of the spectral diameter
            SPDIAM = ISRGHT - ISLEFT;
         } else {
            // Case of bisection
            // Find approximations to the wanted extremal eigenvalues
            ISLEFT = max(GL, W(WBEGIN) - WERR(WBEGIN) - HNDRD * EPS*ABS(W(WBEGIN)- WERR(WBEGIN) ))             ISRGHT = min(GU,W(WEND) + WERR(WEND) + HNDRD * EPS * ABS(W(WEND)+ WERR(WEND)));
         }


         // Decide whether the base representation for the current block
         // L_JBLK D_JBLK L_JBLK^T = T_JBLK - sigma_JBLK I
         // should be on the left or the right end of the current block.
         // The strategy is to shift to the end which is "more populated"
         // Furthermore, decide whether to use DQDS for the computation of
         // the eigenvalue approximations at the end of SLARRE or bisection.
         // dqds is chosen if all eigenvalues are desired or the number of
         // eigenvalues to be computed is large compared to the blocksize.
         if ( ( IRANGE == ALLRNG ) && ( !FORCEB) ) {
            // If all the eigenvalues have to be computed, we use dqd
            USEDQD = true;
            // INDL is the local index of the first eigenvalue to compute
            INDL = 1;
            INDU = IN;
            // MB =  number of eigenvalues to compute
            MB = IN;
            WEND = WBEGIN + MB - 1;
            // Define 1/4 and 3/4 points of the spectrum
            S1 = ISLEFT + FOURTH * SPDIAM;
            S2 = ISRGHT - FOURTH * SPDIAM;
         } else {
            // SLARRD has computed IBLOCK and INDEXW for each eigenvalue
            // approximation.
            // choose sigma
            if ( USEDQD ) {
               S1 = ISLEFT + FOURTH * SPDIAM;
               S2 = ISRGHT - FOURTH * SPDIAM;
            } else {
               TMP = min(ISRGHT,VU) -  max(ISLEFT,VL);
               S1 =  max(ISLEFT,VL) + FOURTH * TMP;
               S2 =  min(ISRGHT,VU) - FOURTH * TMP;
            }
         }

         // Compute the negcount at the 1/4 and 3/4 points
         if (MB > 1) {
            slarrc('T', IN, S1, S2, D(IBEGIN), E(IBEGIN), PIVMIN, CNT, CNT1, CNT2, IINFO);
         }

         if (MB == 1) {
            SIGMA = GL;
            SGNDEF = ONE;
         } else if ( CNT1 - INDL >= INDU - CNT2 ) {
            if ( ( IRANGE == ALLRNG ) && ( !FORCEB) ) {
               SIGMA = max(ISLEFT,GL);
            } else if ( USEDQD ) {
               // use Gerschgorin bound as shift to get pos def matrix
               // for dqds
               SIGMA = ISLEFT;
            } else {
               // use approximation of the first desired eigenvalue of the
               // block as shift
               SIGMA = max(ISLEFT,VL);
            }
            SGNDEF = ONE;
         } else {
            if ( ( IRANGE == ALLRNG ) && ( !FORCEB) ) {
               SIGMA = min(ISRGHT,GU);
            } else if ( USEDQD ) {
               // use Gerschgorin bound as shift to get neg def matrix
               // for dqds
               SIGMA = ISRGHT;
            } else {
               // use approximation of the first desired eigenvalue of the
               // block as shift
               SIGMA = min(ISRGHT,VU);
            }
            SGNDEF = -ONE;
         }


         // An initial SIGMA has been chosen that will be used for computing
         // T - SIGMA I = L D L^T
         // Define the increment TAU of the shift in case the initial shift
         // needs to be refined to obtain a factorization with not too much
         // element growth.
         if ( USEDQD ) {
            // The initial SIGMA was to the outer end of the spectrum
            // the matrix is definite and we need not retreat.
            TAU = SPDIAM*EPS*N + TWO*PIVMIN;
            TAU = max( TAU,TWO*EPS*(SIGMA).abs() );
         } else {
            if (MB > 1) {
               CLWDTH = W(WEND) + WERR(WEND) - W(WBEGIN) - WERR(WBEGIN);
               AVGAP = (CLWDTH / REAL(WEND-WBEGIN)).abs();
               if ( SGNDEF == ONE ) {
                  TAU = HALF*max(WGAP(WBEGIN),AVGAP);
                  TAU = max(TAU,WERR(WBEGIN));
               } else {
                  TAU = HALF*max(WGAP(WEND-1),AVGAP);
                  TAU = max(TAU,WERR(WEND));
               }
            } else {
               TAU = WERR(WBEGIN);
            }
         }

         for (IDUM = 1; IDUM <= MAXTRY; IDUM++) { // 80
            // Compute L D L^T factorization of tridiagonal matrix T - sigma I.
            // Store D in WORK(1:IN), L in WORK(IN+1:2*IN), and reciprocals of
            // pivots in WORK(2*IN+1:3*IN)
            DPIVOT = D( IBEGIN ) - SIGMA;
            WORK[1] = DPIVOT;
            DMAX = ( WORK(1) ).abs();
            J = IBEGIN;
            for (I = 1; I <= IN - 1; I++) { // 70
               WORK[2*IN+I] = ONE / WORK( I );
               TMP = E( J )*WORK( 2*IN+I );
               WORK[IN+I] = TMP;
               DPIVOT = ( D( J+1 )-SIGMA ) - TMP*E( J );
               WORK[I+1] = DPIVOT;
               DMAX = max( DMAX, (DPIVOT).abs() );
               J = J + 1;
            } // 70
            // check for element growth
            if ( DMAX > MAXGROWTH*SPDIAM ) {
               NOREP = true;
            } else {
               NOREP = false;
            }
            if ( USEDQD && !NOREP ) {
               // Ensure the definiteness of the representation
               // All entries of D (of L D L^T) must have the same sign
               for (I = 1; I <= IN; I++) { // 71
                  TMP = SGNDEF*WORK( I );
                  if (TMP < ZERO) NOREP = true ;
               } // 71
            }
            if (NOREP) {
               // Note that in the case of IRANGE=ALLRNG, we use the Gerschgorin
               // shift which makes the matrix definite. So we should end up
               // here really only in the case of IRANGE = VALRNG or INDRNG.
               if ( IDUM == MAXTRY-1 ) {
                  if ( SGNDEF == ONE ) {
                     // The fudged Gerschgorin shift should succeed
                     SIGMA = GL - FUDGE*SPDIAM*EPS*N - FUDGE*TWO*PIVMIN;
                  } else {
                     SIGMA = GU + FUDGE*SPDIAM*EPS*N + FUDGE*TWO*PIVMIN;
                  }
               } else {
                  SIGMA = SIGMA - SGNDEF * TAU;
                  TAU = TWO * TAU;
               }
            } else {
               // an initial RRR is found
               GO TO 83;
            }
         } // 80
         // if the program reaches this point, no base representation could be
         // found in MAXTRY iterations.
         INFO = 2;
         return;

         } // 83
         // At this point, we have found an initial base representation
         // T - SIGMA I = L D L^T with not too much element growth.
         // Store the shift.
         E[IEND] = SIGMA;
         // Store D and L.
         scopy(IN, WORK, 1, D( IBEGIN ), 1 );
         scopy(IN-1, WORK( IN+1 ), 1, E( IBEGIN ), 1 );


         if (MB > 1 ) {

            // Perturb each entry of the base representation by a small
            // (but random) relative amount to overcome difficulties with
            // glued matrices.

            for (I = 1; I <= 4; I++) { // 122
               ISEED[I] = 1;
            } // 122

            slarnv(2, ISEED, 2*IN-1, WORK(1));
            for (I = 1; I <= IN-1; I++) { // 125
               D[IBEGIN+I-1] = D(IBEGIN+I-1)*(ONE+EPS*PERT*WORK(I));
               E[IBEGIN+I-1] = E(IBEGIN+I-1)*(ONE+EPS*PERT*WORK(IN+I));
            } // 125
            D[IEND] = D(IEND)*(ONE+EPS*FOUR*WORK(IN));

         }

         // Don't update the Gerschgorin intervals because keeping track
         // of the updates would be too much work in SLARRV.
         // We update W instead and use it to locate the proper Gerschgorin
         // intervals.

         // Compute the required eigenvalues of L D L' by bisection or dqds
         if ( !USEDQD ) {
            // If SLARRD has been used, shift the eigenvalue approximations
            // according to their representation. This is necessary for
            // a uniform SLARRV since dqds computes eigenvalues of the
            // shifted representation. In SLARRV, W will always hold the
            // UNshifted eigenvalue approximation.
            for (J = WBEGIN; J <= WEND; J++) { // 134
               W[J] = W(J) - SIGMA;
               WERR[J] = WERR(J) + (W(J)).abs() * EPS;
            } // 134
            // call SLARRB to reduce eigenvalue error of the approximations
            // from SLARRD
            for (I = IBEGIN; I <= IEND-1; I++) { // 135
               WORK[I] = D( I ) * E( I )**2;
            } // 135
            // use bisection to find EV from INDL to INDU
            slarrb(IN, D(IBEGIN), WORK(IBEGIN), INDL, INDU, RTOL1, RTOL2, INDL-1, W(WBEGIN), WGAP(WBEGIN), WERR(WBEGIN), WORK( 2*N+1 ), IWORK, PIVMIN, SPDIAM, IN, IINFO );
            if ( IINFO != 0 ) {
               INFO = -4;
               return;
            }
            // SLARRB computes all gaps correctly except for the last one
            // Record distance to VU/GU
            WGAP[WEND] = max( ZERO, ( VU-SIGMA ) - ( W( WEND ) + WERR( WEND ) ) );
            for (I = INDL; I <= INDU; I++) { // 138
               M = M + 1;
               IBLOCK[M] = JBLK;
               INDEXW[M] = I;
            } // 138
         } else {
            // Call dqds to get all eigs (and then possibly delete unwanted
            // eigenvalues).
            // Note that dqds finds the eigenvalues of the L D L^T representation
            // of T to high relative accuracy. High relative accuracy
            // might be lost when the shift of the RRR is subtracted to obtain
            // the eigenvalues of T. However, T is not guaranteed to define its
            // eigenvalues to high relative accuracy anyway.
            // Set RTOL to the order of the tolerance used in SLASQ2
            // This is an ESTIMATED error, the worst case bound is 4*N*EPS
            // which is usually too large and requires unnecessary work to be
            // done by bisection when computing the eigenvectors
            RTOL = LOG(double(IN)) * FOUR * EPS;
            J = IBEGIN;
            for (I = 1; I <= IN - 1; I++) { // 140
               WORK[2*I-1] = ( D( J ) ).abs();
               WORK[2*I] = E( J )*E( J )*WORK( 2*I-1 );
               J = J + 1;
            } // 140
            WORK[2*IN-1] = ( D( IEND ) ).abs();
            WORK[2*IN] = ZERO;
            slasq2(IN, WORK, IINFO );
            if ( IINFO != 0 ) {
               // If IINFO = -5 then an index is part of a tight cluster
               // and should be changed. The index is in IWORK(1) and the
               // gap is in WORK(N+1)
               INFO = -5;
               return;
            } else {
               // Test that all eigenvalues are positive as expected
               for (I = 1; I <= IN; I++) { // 149
                  if ( WORK( I ) < ZERO ) {
                     INFO = -6;
                     return;
                  }
               } // 149
            }
            if ( SGNDEF > ZERO ) {
               for (I = INDL; I <= INDU; I++) { // 150
                  M = M + 1;
                  W[M] = WORK( IN-I+1 );
                  IBLOCK[M] = JBLK;
                  INDEXW[M] = I;
               } // 150
            } else {
               for (I = INDL; I <= INDU; I++) { // 160
                  M = M + 1;
                  W[M] = -WORK( I );
                  IBLOCK[M] = JBLK;
                  INDEXW[M] = I;
               } // 160
            }

            for (I = M - MB + 1; I <= M; I++) { // 165
               // the value of RTOL below should be the tolerance in SLASQ2
               WERR[I] = RTOL * ( W(I) ).abs();
            } // 165
            for (I = M - MB + 1; I <= M - 1; I++) { // 166
               // compute the right gap between the intervals
               WGAP[I] = max( ZERO, W(I+1)-WERR(I+1) - (W(I)+WERR(I)) );
            } // 166
            WGAP[M] = max( ZERO, ( VU-SIGMA ) - ( W( M ) + WERR( M ) ) );
         }
         // proceed with next block
         IBEGIN = IEND + 1;
         WBEGIN = WEND + 1;
      } // 170


      }
