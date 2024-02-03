      void dlaebz(IJOB, NITMAX, N, MMAX, MINP, NBMIN, ABSTOL, RELTOL, PIVMIN, D, E, E2, NVAL, AB, C, MOUT, NAB, WORK, IWORK, INFO ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IJOB, INFO, MINP, MMAX, MOUT, N, NBMIN, NITMAX;
      double             ABSTOL, PIVMIN, RELTOL;
      // ..
      // .. Array Arguments ..
      int                IWORK( * ), NAB( MMAX, * ), NVAL( * );
      double             AB( MMAX, * ), C( * ), D( * ), E( * ), E2( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, TWO, HALF;
      const              ZERO = 0.0, TWO = 2.0, HALF = 1.0 / TWO ;
      // ..
      // .. Local Scalars ..
      int                ITMP1, ITMP2, J, JI, JIT, JP, KF, KFNEW, KL, KLNEW;
      double             TMP1, TMP2;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Check for Errors

      INFO = 0;
      if ( IJOB < 1 || IJOB > 3 ) {
         INFO = -1;
         return;
      }

      // Initialize NAB

      if ( IJOB == 1 ) {

         // Compute the number of eigenvalues in the initial intervals.

         MOUT = 0;
         for (JI = 1; JI <= MINP; JI++) { // 30
            for (JP = 1; JP <= 2; JP++) { // 20
               TMP1 = D( 1 ) - AB( JI, JP );
               if( ABS( TMP1 ) < PIVMIN ) TMP1 = -PIVMIN;
               NAB( JI, JP ) = 0;
               if (TMP1 <= ZERO) NAB( JI, JP ) = 1;

               for (J = 2; J <= N; J++) { // 10
                  TMP1 = D( J ) - E2( J-1 ) / TMP1 - AB( JI, JP );
                  if( ABS( TMP1 ) < PIVMIN ) TMP1 = -PIVMIN;
                  IF( TMP1 <= ZERO ) NAB( JI, JP ) = NAB( JI, JP ) + 1;
               } // 10
            } // 20
            MOUT = MOUT + NAB( JI, 2 ) - NAB( JI, 1 );
         } // 30
         return;
      }

      // Initialize for loop

      // KF and KL have the following meaning:
         // Intervals 1,...,KF-1 have converged.
         // Intervals KF,...,KL  still need to be refined.

      KF = 1;
      KL = MINP;

      // If IJOB=2, initialize C.
      // If IJOB=3, use the user-supplied starting point.

      if ( IJOB == 2 ) {
         for (JI = 1; JI <= MINP; JI++) { // 40
            C( JI ) = HALF*( AB( JI, 1 )+AB( JI, 2 ) );
         } // 40
      }

      // Iteration loop

      for (JIT = 1; JIT <= NITMAX; JIT++) { // 130

         // Loop over intervals

         if ( KL-KF+1 >= NBMIN && NBMIN > 0 ) {

            // Begin of Parallel Version of the loop

            for (JI = KF; JI <= KL; JI++) { // 60

               // Compute N(c), the number of eigenvalues less than c

               WORK( JI ) = D( 1 ) - C( JI );
               IWORK( JI ) = 0;
               if ( WORK( JI ) <= PIVMIN ) {
                  IWORK( JI ) = 1;
                  WORK( JI ) = min( WORK( JI ), -PIVMIN );
               }

               for (J = 2; J <= N; J++) { // 50
                  WORK( JI ) = D( J ) - E2( J-1 ) / WORK( JI ) - C( JI );
                  if ( WORK( JI ) <= PIVMIN ) {
                     IWORK( JI ) = IWORK( JI ) + 1;
                     WORK( JI ) = min( WORK( JI ), -PIVMIN );
                  }
               } // 50
            } // 60

            if ( IJOB <= 2 ) {

               // IJOB=2: Choose all intervals containing eigenvalues.

               KLNEW = KL;
               for (JI = KF; JI <= KL; JI++) { // 70

                  // Insure that N(w) is monotone

                  IWORK( JI ) = min( NAB( JI, 2 ), max( NAB( JI, 1 ), IWORK( JI ) ) );

                  // Update the Queue -- add intervals if both halves
                  // contain eigenvalues.

                  if ( IWORK( JI ) == NAB( JI, 2 ) ) {

                     // No eigenvalue in the upper interval:
                     // just use the lower interval.

                     AB( JI, 2 ) = C( JI );

                  } else if ( IWORK( JI ) == NAB( JI, 1 ) ) {

                     // No eigenvalue in the lower interval:
                     // just use the upper interval.

                     AB( JI, 1 ) = C( JI );
                  } else {
                     KLNEW = KLNEW + 1;
                     if ( KLNEW <= MMAX ) {

                        // Eigenvalue in both intervals -- add upper to
                        // queue.

                        AB( KLNEW, 2 ) = AB( JI, 2 );
                        NAB( KLNEW, 2 ) = NAB( JI, 2 );
                        AB( KLNEW, 1 ) = C( JI );
                        NAB( KLNEW, 1 ) = IWORK( JI );
                        AB( JI, 2 ) = C( JI );
                        NAB( JI, 2 ) = IWORK( JI );
                     } else {
                        INFO = MMAX + 1;
                     }
                  }
               } // 70
               if (INFO != 0) return;
               KL = KLNEW;
            } else {

               // IJOB=3: Binary search.  Keep only the interval containing
                       // w   s.t. N(w) = NVAL

               for (JI = KF; JI <= KL; JI++) { // 80
                  if ( IWORK( JI ) <= NVAL( JI ) ) {
                     AB( JI, 1 ) = C( JI );
                     NAB( JI, 1 ) = IWORK( JI );
                  }
                  if ( IWORK( JI ) >= NVAL( JI ) ) {
                     AB( JI, 2 ) = C( JI );
                     NAB( JI, 2 ) = IWORK( JI );
                  }
               } // 80
            }

         } else {

            // End of Parallel Version of the loop

            // Begin of Serial Version of the loop

            KLNEW = KL;
            for (JI = KF; JI <= KL; JI++) { // 100

               // Compute N(w), the number of eigenvalues less than w

               TMP1 = C( JI );
               TMP2 = D( 1 ) - TMP1;
               ITMP1 = 0;
               if ( TMP2 <= PIVMIN ) {
                  ITMP1 = 1;
                  TMP2 = min( TMP2, -PIVMIN );
               }

               for (J = 2; J <= N; J++) { // 90
                  TMP2 = D( J ) - E2( J-1 ) / TMP2 - TMP1;
                  if ( TMP2 <= PIVMIN ) {
                     ITMP1 = ITMP1 + 1;
                     TMP2 = min( TMP2, -PIVMIN );
                  }
               } // 90

               if ( IJOB <= 2 ) {

                  // IJOB=2: Choose all intervals containing eigenvalues.

                  // Insure that N(w) is monotone

                  ITMP1 = min( NAB( JI, 2 ), max( NAB( JI, 1 ), ITMP1 ) );

                  // Update the Queue -- add intervals if both halves
                  // contain eigenvalues.

                  if ( ITMP1 == NAB( JI, 2 ) ) {

                     // No eigenvalue in the upper interval:
                     // just use the lower interval.

                     AB( JI, 2 ) = TMP1;

                  } else if ( ITMP1 == NAB( JI, 1 ) ) {

                     // No eigenvalue in the lower interval:
                     // just use the upper interval.

                     AB( JI, 1 ) = TMP1;
                  } else if ( KLNEW < MMAX ) {

                     // Eigenvalue in both intervals -- add upper to queue.

                     KLNEW = KLNEW + 1;
                     AB( KLNEW, 2 ) = AB( JI, 2 );
                     NAB( KLNEW, 2 ) = NAB( JI, 2 );
                     AB( KLNEW, 1 ) = TMP1;
                     NAB( KLNEW, 1 ) = ITMP1;
                     AB( JI, 2 ) = TMP1;
                     NAB( JI, 2 ) = ITMP1;
                  } else {
                     INFO = MMAX + 1;
                     return;
                  }
               } else {

                  // IJOB=3: Binary search.  Keep only the interval
                          // containing  w  s.t. N(w) = NVAL

                  if ( ITMP1 <= NVAL( JI ) ) {
                     AB( JI, 1 ) = TMP1;
                     NAB( JI, 1 ) = ITMP1;
                  }
                  if ( ITMP1 >= NVAL( JI ) ) {
                     AB( JI, 2 ) = TMP1;
                     NAB( JI, 2 ) = ITMP1;
                  }
               }
            } // 100
            KL = KLNEW;

         }

         // Check for convergence

         KFNEW = KF;
         for (JI = KF; JI <= KL; JI++) { // 110
            TMP1 = ABS( AB( JI, 2 )-AB( JI, 1 ) );
            TMP2 = max( ABS( AB( JI, 2 ) ), ABS( AB( JI, 1 ) ) );
            if ( TMP1 < max( ABSTOL, PIVMIN, RELTOL*TMP2 ) || NAB( JI, 1 ) >= NAB( JI, 2 ) ) {

               // Converged -- Swap with position KFNEW,
                            // then increment KFNEW

               if ( JI > KFNEW ) {
                  TMP1 = AB( JI, 1 );
                  TMP2 = AB( JI, 2 );
                  ITMP1 = NAB( JI, 1 );
                  ITMP2 = NAB( JI, 2 );
                  AB( JI, 1 ) = AB( KFNEW, 1 );
                  AB( JI, 2 ) = AB( KFNEW, 2 );
                  NAB( JI, 1 ) = NAB( KFNEW, 1 );
                  NAB( JI, 2 ) = NAB( KFNEW, 2 );
                  AB( KFNEW, 1 ) = TMP1;
                  AB( KFNEW, 2 ) = TMP2;
                  NAB( KFNEW, 1 ) = ITMP1;
                  NAB( KFNEW, 2 ) = ITMP2;
                  if ( IJOB == 3 ) {
                     ITMP1 = NVAL( JI );
                     NVAL( JI ) = NVAL( KFNEW );
                     NVAL( KFNEW ) = ITMP1;
                  }
               }
               KFNEW = KFNEW + 1;
            }
         } // 110
         KF = KFNEW;

         // Choose Midpoints

         for (JI = KF; JI <= KL; JI++) { // 120
            C( JI ) = HALF*( AB( JI, 1 )+AB( JI, 2 ) );
         } // 120

         // If no more intervals to refine, quit.

         if (KF > KL) GO TO 140;
      } // 130

      // Converged

      } // 140
      INFO = max( KL+1-KF, 0 );
      MOUT = KL;

      return;

      // End of DLAEBZ

      }
