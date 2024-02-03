      SUBROUTINE SLAEBZ( IJOB, NITMAX, N, MMAX, MINP, NBMIN, ABSTOL, RELTOL, PIVMIN, D, E, E2, NVAL, AB, C, MOUT, NAB, WORK, IWORK, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IJOB, INFO, MINP, MMAX, MOUT, N, NBMIN, NITMAX;
      REAL               ABSTOL, PIVMIN, RELTOL
      // ..
      // .. Array Arguments ..
      int                IWORK( * ), NAB( MMAX, * ), NVAL( * );
      REAL               AB( MMAX, * ), C( * ), D( * ), E( * ), E2( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, TWO, HALF
      const              ZERO = 0.0E0, TWO = 2.0E0, HALF = 1.0E0 / TWO ;
      // ..
      // .. Local Scalars ..
      int                ITMP1, ITMP2, J, JI, JIT, JP, KF, KFNEW, KL, KLNEW;
      REAL               TMP1, TMP2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Check for Errors

      INFO = 0
      IF( IJOB.LT.1 .OR. IJOB.GT.3 ) THEN
         INFO = -1
         RETURN
      END IF

      // Initialize NAB

      IF( IJOB.EQ.1 ) THEN

         // Compute the number of eigenvalues in the initial intervals.

         MOUT = 0
         DO 30 JI = 1, MINP
            DO 20 JP = 1, 2
               TMP1 = D( 1 ) - AB( JI, JP )
               IF( ABS( TMP1 ).LT.PIVMIN ) TMP1 = -PIVMIN
               NAB( JI, JP ) = 0
               IF( TMP1.LE.ZERO ) NAB( JI, JP ) = 1

               DO 10 J = 2, N
                  TMP1 = D( J ) - E2( J-1 ) / TMP1 - AB( JI, JP )
                  IF( ABS( TMP1 ).LT.PIVMIN ) TMP1 = -PIVMIN                   IF( TMP1.LE.ZERO ) NAB( JI, JP ) = NAB( JI, JP ) + 1
   10          CONTINUE
   20       CONTINUE
            MOUT = MOUT + NAB( JI, 2 ) - NAB( JI, 1 )
   30    CONTINUE
         RETURN
      END IF

      // Initialize for loop

      // KF and KL have the following meaning:
         // Intervals 1,...,KF-1 have converged.
         // Intervals KF,...,KL  still need to be refined.

      KF = 1
      KL = MINP

      // If IJOB=2, initialize C.
      // If IJOB=3, use the user-supplied starting point.

      IF( IJOB.EQ.2 ) THEN
         DO 40 JI = 1, MINP
            C( JI ) = HALF*( AB( JI, 1 )+AB( JI, 2 ) )
   40    CONTINUE
      END IF

      // Iteration loop

      DO 130 JIT = 1, NITMAX

         // Loop over intervals

         IF( KL-KF+1.GE.NBMIN .AND. NBMIN.GT.0 ) THEN

            // Begin of Parallel Version of the loop

            DO 60 JI = KF, KL

               // Compute N(c), the number of eigenvalues less than c

               WORK( JI ) = D( 1 ) - C( JI )
               IWORK( JI ) = 0
               IF( WORK( JI ).LE.PIVMIN ) THEN
                  IWORK( JI ) = 1
                  WORK( JI ) = MIN( WORK( JI ), -PIVMIN )
               END IF

               DO 50 J = 2, N
                  WORK( JI ) = D( J ) - E2( J-1 ) / WORK( JI ) - C( JI )
                  IF( WORK( JI ).LE.PIVMIN ) THEN
                     IWORK( JI ) = IWORK( JI ) + 1
                     WORK( JI ) = MIN( WORK( JI ), -PIVMIN )
                  END IF
   50          CONTINUE
   60       CONTINUE

            IF( IJOB.LE.2 ) THEN

               // IJOB=2: Choose all intervals containing eigenvalues.

               KLNEW = KL
               DO 70 JI = KF, KL

                  // Insure that N(w) is monotone

                  IWORK( JI ) = MIN( NAB( JI, 2 ), MAX( NAB( JI, 1 ), IWORK( JI ) ) )

                  // Update the Queue -- add intervals if both halves
                  // contain eigenvalues.

                  IF( IWORK( JI ).EQ.NAB( JI, 2 ) ) THEN

                     // No eigenvalue in the upper interval:
                     // just use the lower interval.

                     AB( JI, 2 ) = C( JI )

                  ELSE IF( IWORK( JI ).EQ.NAB( JI, 1 ) ) THEN

                     // No eigenvalue in the lower interval:
                     // just use the upper interval.

                     AB( JI, 1 ) = C( JI )
                  } else {
                     KLNEW = KLNEW + 1
                     IF( KLNEW.LE.MMAX ) THEN

                        // Eigenvalue in both intervals -- add upper to
                        // queue.

                        AB( KLNEW, 2 ) = AB( JI, 2 )
                        NAB( KLNEW, 2 ) = NAB( JI, 2 )
                        AB( KLNEW, 1 ) = C( JI )
                        NAB( KLNEW, 1 ) = IWORK( JI )
                        AB( JI, 2 ) = C( JI )
                        NAB( JI, 2 ) = IWORK( JI )
                     } else {
                        INFO = MMAX + 1
                     END IF
                  END IF
   70          CONTINUE
               IF( INFO.NE.0 ) RETURN
               KL = KLNEW
            } else {

               // IJOB=3: Binary search.  Keep only the interval containing
                       // w   s.t. N(w) = NVAL

               DO 80 JI = KF, KL
                  IF( IWORK( JI ).LE.NVAL( JI ) ) THEN
                     AB( JI, 1 ) = C( JI )
                     NAB( JI, 1 ) = IWORK( JI )
                  END IF
                  IF( IWORK( JI ).GE.NVAL( JI ) ) THEN
                     AB( JI, 2 ) = C( JI )
                     NAB( JI, 2 ) = IWORK( JI )
                  END IF
   80          CONTINUE
            END IF

         } else {

            // End of Parallel Version of the loop

            // Begin of Serial Version of the loop

            KLNEW = KL
            DO 100 JI = KF, KL

               // Compute N(w), the number of eigenvalues less than w

               TMP1 = C( JI )
               TMP2 = D( 1 ) - TMP1
               ITMP1 = 0
               IF( TMP2.LE.PIVMIN ) THEN
                  ITMP1 = 1
                  TMP2 = MIN( TMP2, -PIVMIN )
               END IF

               DO 90 J = 2, N
                  TMP2 = D( J ) - E2( J-1 ) / TMP2 - TMP1
                  IF( TMP2.LE.PIVMIN ) THEN
                     ITMP1 = ITMP1 + 1
                     TMP2 = MIN( TMP2, -PIVMIN )
                  END IF
   90          CONTINUE

               IF( IJOB.LE.2 ) THEN

                  // IJOB=2: Choose all intervals containing eigenvalues.

                  // Insure that N(w) is monotone

                  ITMP1 = MIN( NAB( JI, 2 ), MAX( NAB( JI, 1 ), ITMP1 ) )

                  // Update the Queue -- add intervals if both halves
                  // contain eigenvalues.

                  IF( ITMP1.EQ.NAB( JI, 2 ) ) THEN

                     // No eigenvalue in the upper interval:
                     // just use the lower interval.

                     AB( JI, 2 ) = TMP1

                  ELSE IF( ITMP1.EQ.NAB( JI, 1 ) ) THEN

                     // No eigenvalue in the lower interval:
                     // just use the upper interval.

                     AB( JI, 1 ) = TMP1
                  ELSE IF( KLNEW.LT.MMAX ) THEN

                     // Eigenvalue in both intervals -- add upper to queue.

                     KLNEW = KLNEW + 1
                     AB( KLNEW, 2 ) = AB( JI, 2 )
                     NAB( KLNEW, 2 ) = NAB( JI, 2 )
                     AB( KLNEW, 1 ) = TMP1
                     NAB( KLNEW, 1 ) = ITMP1
                     AB( JI, 2 ) = TMP1
                     NAB( JI, 2 ) = ITMP1
                  } else {
                     INFO = MMAX + 1
                     RETURN
                  END IF
               } else {

                  // IJOB=3: Binary search.  Keep only the interval
                          // containing  w  s.t. N(w) = NVAL

                  IF( ITMP1.LE.NVAL( JI ) ) THEN
                     AB( JI, 1 ) = TMP1
                     NAB( JI, 1 ) = ITMP1
                  END IF
                  IF( ITMP1.GE.NVAL( JI ) ) THEN
                     AB( JI, 2 ) = TMP1
                     NAB( JI, 2 ) = ITMP1
                  END IF
               END IF
  100       CONTINUE
            KL = KLNEW

         END IF

         // Check for convergence

         KFNEW = KF
         DO 110 JI = KF, KL
            TMP1 = ABS( AB( JI, 2 )-AB( JI, 1 ) )
            TMP2 = MAX( ABS( AB( JI, 2 ) ), ABS( AB( JI, 1 ) ) )
            IF( TMP1.LT.MAX( ABSTOL, PIVMIN, RELTOL*TMP2 ) .OR. NAB( JI, 1 ).GE.NAB( JI, 2 ) ) THEN

               // Converged -- Swap with position KFNEW,
                           t // hen increment KFNEW

               IF( JI.GT.KFNEW ) THEN
                  TMP1 = AB( JI, 1 )
                  TMP2 = AB( JI, 2 )
                  ITMP1 = NAB( JI, 1 )
                  ITMP2 = NAB( JI, 2 )
                  AB( JI, 1 ) = AB( KFNEW, 1 )
                  AB( JI, 2 ) = AB( KFNEW, 2 )
                  NAB( JI, 1 ) = NAB( KFNEW, 1 )
                  NAB( JI, 2 ) = NAB( KFNEW, 2 )
                  AB( KFNEW, 1 ) = TMP1
                  AB( KFNEW, 2 ) = TMP2
                  NAB( KFNEW, 1 ) = ITMP1
                  NAB( KFNEW, 2 ) = ITMP2
                  IF( IJOB.EQ.3 ) THEN
                     ITMP1 = NVAL( JI )
                     NVAL( JI ) = NVAL( KFNEW )
                     NVAL( KFNEW ) = ITMP1
                  END IF
               END IF
               KFNEW = KFNEW + 1
            END IF
  110    CONTINUE
         KF = KFNEW

         // Choose Midpoints

         DO 120 JI = KF, KL
            C( JI ) = HALF*( AB( JI, 1 )+AB( JI, 2 ) )
  120    CONTINUE

         // If no more intervals to refine, quit.

         IF( KF.GT.KL ) GO TO 140
  130 CONTINUE

      // Converged

  140 CONTINUE
      INFO = MAX( KL+1-KF, 0 )
      MOUT = KL

      RETURN

      // End of SLAEBZ

      }
