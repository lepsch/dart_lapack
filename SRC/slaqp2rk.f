      SUBROUTINE SLAQP2RK( M, N, NRHS, IOFFSET, KMAX, ABSTOL, RELTOL, KP1, MAXC2NRM, A, LDA, K, MAXC2NRMK, RELMAXC2NRMK, JPIV, TAU, VN1, VN2, WORK, INFO )
      IMPLICIT NONE

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, IOFFSET, KP1, K, KMAX, LDA, M, N, NRHS;
      REAL               ABSTOL, MAXC2NRM, MAXC2NRMK, RELMAXC2NRMK, RELTOL
      // ..
      // .. Array Arguments ..
      int                JPIV( * );
      REAL               A( LDA, * ), TAU( * ), VN1( * ), VN2( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                I, ITEMP, J, JMAXC2NRM, KK, KP, MINMNFACT, MINMNUPDT;
      REAL               AIKK, HUGEVAL, TEMP, TEMP2, TOL3Z
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLARF, SLARFG, SSWAP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, SQRT
      // ..
      // .. External Functions ..
      bool               SISNAN;
      int                ISAMAX;
      REAL               SLAMCH, SNRM2
      // EXTERNAL SISNAN, SLAMCH, ISAMAX, SNRM2
      // ..
      // .. Executable Statements ..

      // Initialize INFO

      INFO = 0

      // MINMNFACT in the smallest dimension of the submatrix
      // A(IOFFSET+1:M,1:N) to be factorized.

      // MINMNUPDT is the smallest dimension
      // of the subarray A(IOFFSET+1:M,1:N+NRHS) to be udated, which
      // contains the submatrices A(IOFFSET+1:M,1:N) and
      // B(IOFFSET+1:M,1:NRHS) as column blocks.

      MINMNFACT = MIN( M-IOFFSET, N )
      MINMNUPDT = MIN( M-IOFFSET, N+NRHS )
      KMAX = MIN( KMAX, MINMNFACT )
      TOL3Z = SQRT( SLAMCH( 'Epsilon' ) )
      HUGEVAL = SLAMCH( 'Overflow' )

      // Compute the factorization, KK is the lomn loop index.

      for (KK = 1; KK <= KMAX; KK++) {

         I = IOFFSET + KK

         if ( I.EQ.1 ) {

            // ============================================================

            // We are at the first column of the original whole matrix A,
            // therefore we use the computed KP1 and MAXC2NRM from the
            // main routine.


            KP = KP1

            // ============================================================

         } else {

            // ============================================================

            // Determine the pivot column in KK-th step, i.e. the index
            // of the column with the maximum 2-norm in the
            // submatrix A(I:M,K:N).

            KP = ( KK-1 ) + ISAMAX( N-KK+1, VN1( KK ), 1 )

            // Determine the maximum column 2-norm and the relative maximum
            // column 2-norm of the submatrix A(I:M,KK:N) in step KK.
            // RELMAXC2NRMK  will be computed later, after somecondition
            // checks on MAXC2NRMK.

            MAXC2NRMK = VN1( KP )

            // ============================================================

            // Check if the submatrix A(I:M,KK:N) contains NaN, and set
            // INFO parameter to the column number, where the first NaN
            // is found and return from the routine.
            // We need to check the condition only if the
            // column index (same as row index) of the original whole
            // matrix is larger than 1, since the condition for whole
            // original matrix is checked in the main routine.

            if ( SISNAN( MAXC2NRMK ) ) {

               // Set K, the number of factorized columns.
               // that are not zero.

                K = KK - 1
                INFO = K + KP

                // Set RELMAXC2NRMK to NaN.

                RELMAXC2NRMK = MAXC2NRMK

                // Array TAU(K+1:MINMNFACT) is not set and contains
                // undefined elements.

               RETURN
            }

            // ============================================================

            // Quick return, if the submatrix A(I:M,KK:N) is
            // a zero matrix.
            // We need to check the condition only if the
            // column index (same as row index) of the original whole
            // matrix is larger than 1, since the condition for whole
            // original matrix is checked in the main routine.

            if ( MAXC2NRMK.EQ.ZERO ) {

               // Set K, the number of factorized columns.
               // that are not zero.

               K = KK - 1
               RELMAXC2NRMK = ZERO

               // Set TAUs corresponding to the columns that were not
               // factorized to ZERO, i.e. set TAU(KK:MINMNFACT) to ZERO.

               for (J = KK; J <= MINMNFACT; J++) {
                  TAU( J ) = ZERO
               }

               // Return from the routine.

               RETURN

            }

            // ============================================================

            // Check if the submatrix A(I:M,KK:N) contains Inf,
            // set INFO parameter to the column number, where
            // the first Inf is found plus N, and continue
            // the computation.
            // We need to check the condition only if the
            // column index (same as row index) of the original whole
            // matrix is larger than 1, since the condition for whole
            // original matrix is checked in the main routine.

            if ( INFO.EQ.0 .AND. MAXC2NRMK.GT.HUGEVAL ) {
               INFO = N + KK - 1 + KP
            }

            // ============================================================

            // Test for the second and third stopping criteria.
            // NOTE: There is no need to test for ABSTOL >= ZERO, since
            // MAXC2NRMK is non-negative. Similarly, there is no need
            // to test for RELTOL >= ZERO, since RELMAXC2NRMK is
            // non-negative.
            // We need to check the condition only if the
            // column index (same as row index) of the original whole
            // matrix is larger than 1, since the condition for whole
            // original matrix is checked in the main routine.

            RELMAXC2NRMK =  MAXC2NRMK / MAXC2NRM

            if ( MAXC2NRMK.LE.ABSTOL .OR. RELMAXC2NRMK.LE.RELTOL ) {

               // Set K, the number of factorized columns.

               K = KK - 1

               // Set TAUs corresponding to the columns that were not
               // factorized to ZERO, i.e. set TAU(KK:MINMNFACT) to ZERO.

               for (J = KK; J <= MINMNFACT; J++) {
                  TAU( J ) = ZERO
               }

               // Return from the routine.

               RETURN

            }

            // ============================================================

            // End ELSE of IF(I.EQ.1)

         }

         // ===============================================================

         // If the pivot column is not the first column of the
         // subblock A(1:M,KK:N):
         // 1) swap the KK-th column and the KP-th pivot column
            // in A(1:M,1:N);
         // 2) copy the KK-th element into the KP-th element of the partial
            // and exact 2-norm vectors VN1 and VN2. ( Swap is not needed
            // for VN1 and VN2 since we use the element with the index
            // larger than KK in the next loop step.)
         // 3) Save the pivot interchange with the indices relative to the
            // the original matrix A, not the block A(1:M,1:N).

         if ( KP.NE.KK ) {
            sswap(M, A( 1, KP ), 1, A( 1, KK ), 1 );
            VN1( KP ) = VN1( KK )
            VN2( KP ) = VN2( KK )
            ITEMP = JPIV( KP )
            JPIV( KP ) = JPIV( KK )
            JPIV( KK ) = ITEMP
         }

         // Generate elementary reflector H(KK) using the column A(I:M,KK),
         // if the column has more than one element, otherwise
         // the elementary reflector would be an identity matrix,
         // and TAU(KK) = ZERO.

         if ( I.LT.M ) {
            slarfg(M-I+1, A( I, KK ), A( I+1, KK ), 1, TAU( KK ) );
         } else {
            TAU( KK ) = ZERO
         }

         // Check if TAU(KK) contains NaN, set INFO parameter
         // to the column number where NaN is found and return from
         // the routine.
         // NOTE: There is no need to check TAU(KK) for Inf,
         // since SLARFG cannot produce TAU(KK) or Householder vector
         // below the diagonal containing Inf. Only BETA on the diagonal,
         // returned by SLARFG can contain Inf, which requires
         // TAU(KK) to contain NaN. Therefore, this case of generating Inf
         // by SLARFG is covered by checking TAU(KK) for NaN.

         if ( SISNAN( TAU(KK) ) ) {
            K = KK - 1
            INFO = KK

            // Set MAXC2NRMK and  RELMAXC2NRMK to NaN.

            MAXC2NRMK = TAU( KK )
            RELMAXC2NRMK = TAU( KK )

            // Array TAU(KK:MINMNFACT) is not set and contains
            // undefined elements, except the first element TAU(KK) = NaN.

            RETURN
         }

         // Apply H(KK)**T to A(I:M,KK+1:N+NRHS) from the left.
         // ( If M >= N, then at KK = N there is no residual matrix,
          // i.e. no columns of A to update, only columns of B.
          // If M < N, then at KK = M-IOFFSET, I = M and we have a
          // one-row residual matrix in A and the elementary
          // reflector is a unit matrix, TAU(KK) = ZERO, i.e. no update
          // is needed for the residual matrix in A and the
          // right-hand-side-matrix in B.
          // Therefore, we update only if
          // KK < MINMNUPDT = min(M-IOFFSET, N+NRHS)
          // condition is satisfied, not only KK < N+NRHS )

         if ( KK.LT.MINMNUPDT ) {
            AIKK = A( I, KK )
            A( I, KK ) = ONE
            slarf('Left', M-I+1, N+NRHS-KK, A( I, KK ), 1, TAU( KK ), A( I, KK+1 ), LDA, WORK( 1 ) );
            A( I, KK ) = AIKK
         }

         if ( KK.LT.MINMNFACT ) {

            // Update the partial column 2-norms for the residual matrix,
            // only if the residual matrix A(I+1:M,KK+1:N) exists, i.e.
            // when KK < min(M-IOFFSET, N).

            for (J = KK + 1; J <= N; J++) {
               if ( VN1( J ).NE.ZERO ) {

                  // NOTE: The following lines follow from the analysis in
                  // Lapack Working Note 176.

                  TEMP = ONE - ( ABS( A( I, J ) ) / VN1( J ) )**2
                  TEMP = MAX( TEMP, ZERO )
                  TEMP2 = TEMP*( VN1( J ) / VN2( J ) )**2
                  if ( TEMP2 .LE. TOL3Z ) {

                     // Compute the column 2-norm for the partial
                     // column A(I+1:M,J) by explicitly computing it,
                     // and store it in both partial 2-norm vector VN1
                     // and exact column 2-norm vector VN2.

                     VN1( J ) = SNRM2( M-I, A( I+1, J ), 1 )
                     VN2( J ) = VN1( J )

                  } else {

                     // Update the column 2-norm for the partial
                     // column A(I+1:M,J) by removing one
                     // element A(I,J) and store it in partial
                     // 2-norm vector VN1.

                     VN1( J ) = VN1( J )*SQRT( TEMP )

                  }
               }
            }

         }

      // End factorization loop

      }

      // If we reached this point, all colunms have been factorized,
      // i.e. no condition was triggered to exit the routine.
      // Set the number of factorized columns.

      K = KMAX

      // We reached the end of the loop, i.e. all KMAX columns were
      // factorized, we need to set MAXC2NRMK and RELMAXC2NRMK before
      // we return.

      if ( K.LT.MINMNFACT ) {

         JMAXC2NRM = K + ISAMAX( N-K, VN1( K+1 ), 1 )
         MAXC2NRMK = VN1( JMAXC2NRM )

         if ( K.EQ.0 ) {
            RELMAXC2NRMK = ONE
         } else {
            RELMAXC2NRMK = MAXC2NRMK / MAXC2NRM
         }

      } else {
         MAXC2NRMK = ZERO
         RELMAXC2NRMK = ZERO
      }

      // We reached the end of the loop, i.e. all KMAX columns were
      // factorized, set TAUs corresponding to the columns that were
      // not factorized to ZERO, i.e. TAU(K+1:MINMNFACT) set to ZERO.

      for (J = K + 1; J <= MINMNFACT; J++) {
         TAU( J ) = ZERO
      }

      RETURN

      // End of SLAQP2RK

      }
