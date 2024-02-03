      SUBROUTINE CLAQP3RK( M, N, NRHS, IOFFSET, NB, ABSTOL, RELTOL, KP1, MAXC2NRM, A, LDA, DONE, KB, MAXC2NRMK, RELMAXC2NRMK, JPIV, TAU, VN1, VN2, AUXV, F, LDF, IWORK, INFO )
      IMPLICIT NONE

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               DONE;
      int                INFO, IOFFSET, KB, KP1, LDA, LDF, M, N, NB, NRHS;
      REAL               ABSTOL, MAXC2NRM, MAXC2NRMK, RELMAXC2NRMK, RELTOL;
      // ..
      // .. Array Arguments ..
      int                IWORK( * ), JPIV( * );
      REAL               VN1( * ), VN2( * )
      COMPLEX            A( LDA, * ), AUXV( * ), F( LDF, * ), TAU( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      COMPLEX            CZERO, CONE
      const              CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      int                ITEMP, J, K, MINMNFACT, MINMNUPDT, LSTICC, KP, I, IF;
      REAL               HUGEVAL, TAUNAN, TEMP, TEMP2, TOL3Z
      COMPLEX            AIK
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CGEMV, CLARFG, CSWAP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, REAL, CONJG, AIMAG, MAX, MIN, SQRT
      // ..
      // .. External Functions ..
      bool               SISNAN;
      int                ISAMAX;
      REAL               SLAMCH, SCNRM2
      // EXTERNAL SISNAN, SLAMCH, ISAMAX, SCNRM2
      // ..
      // .. Executable Statements ..

      // Initialize INFO

      INFO = 0

      // MINMNFACT in the smallest dimension of the submatrix
      // A(IOFFSET+1:M,1:N) to be factorized.

      MINMNFACT = MIN( M-IOFFSET, N )
      MINMNUPDT = MIN( M-IOFFSET, N+NRHS )
      NB = MIN( NB, MINMNFACT )
      TOL3Z = SQRT( SLAMCH( 'Epsilon' ) )
      HUGEVAL = SLAMCH( 'Overflow' )

      // Compute factorization in a while loop over NB columns,
      // K is the column index in the block A(1:M,1:N).

      K = 0
      LSTICC = 0
      DONE = .FALSE.

      DO WHILE ( K.LT.NB .AND. LSTICC.EQ.0 )
         K = K + 1
         I = IOFFSET + K

         if ( I.EQ.1 ) {

            // We are at the first column of the original whole matrix A_orig,
            // therefore we use the computed KP1 and MAXC2NRM from the
            // main routine.

            KP = KP1

         } else {

            // Determine the pivot column in K-th step, i.e. the index
            // of the column with the maximum 2-norm in the
            // submatrix A(I:M,K:N).

            KP = ( K-1 ) + ISAMAX( N-K+1, VN1( K ), 1 )

            // Determine the maximum column 2-norm and the relative maximum
            // column 2-norm of the submatrix A(I:M,K:N) in step K.

            MAXC2NRMK = VN1( KP )

            // ============================================================

            // Check if the submatrix A(I:M,K:N) contains NaN, set
            // INFO parameter to the column number, where the first NaN
            // is found and return from the routine.
            // We need to check the condition only if the
            // column index (same as row index) of the original whole
            // matrix is larger than 1, since the condition for whole
            // original matrix is checked in the main routine.

            if ( SISNAN( MAXC2NRMK ) ) {

               DONE = .TRUE.

               // Set KB, the number of factorized partial columns
                       // that are non-zero in each step in the block,
                       // i.e. the rank of the factor R.
               // Set IF, the number of processed rows in the block, which
                       // is the same as the number of processed rows in
                       // the original whole matrix A_orig.

               KB = K - 1
               IF = I - 1
               INFO = KB + KP

               // Set RELMAXC2NRMK to NaN.

               RELMAXC2NRMK = MAXC2NRMK

               // There is no need to apply the block reflector to the
               // residual of the matrix A stored in A(KB+1:M,KB+1:N),
               // since the submatrix contains NaN and we stop
               // the computation.
               // But, we need to apply the block reflector to the residual
               // right hand sides stored in A(KB+1:M,N+1:N+NRHS), if the
               // residual right hand sides exist.  This occurs
               // when ( NRHS != 0 AND KB <= (M-IOFFSET) ):

               // A(I+1:M,N+1:N+NRHS) := A(I+1:M,N+1:N+NRHS) -
                                // A(I+1:M,1:KB) * F(N+1:N+NRHS,1:KB)**H.

               if ( NRHS.GT.0 .AND. KB.LT.(M-IOFFSET) ) {
                  cgemm('No transpose', 'Conjugate transpose', M-IF, NRHS, KB, -CONE, A( IF+1, 1 ), LDA, F( N+1, 1 ), LDF, CONE, A( IF+1, N+1 ), LDA );
               }

               // There is no need to recompute the 2-norm of the
               // difficult columns, since we stop the factorization.

               // Array TAU(KF+1:MINMNFACT) is not set and contains
               // undefined elements.

               // Return from the routine.

               RETURN
            }

            // Quick return, if the submatrix A(I:M,K:N) is
            // a zero matrix. We need to check it only if the column index
            // (same as row index) is larger than 1, since the condition
            // for the whole original matrix A_orig is checked in the main
            // routine.

            if ( MAXC2NRMK.EQ.ZERO ) {

               DONE = .TRUE.

               // Set KB, the number of factorized partial columns
                       // that are non-zero in each step in the block,
                       // i.e. the rank of the factor R.
               // Set IF, the number of processed rows in the block, which
                       // is the same as the number of processed rows in
                       // the original whole matrix A_orig.

               KB = K - 1
               IF = I - 1
               RELMAXC2NRMK = ZERO

               // There is no need to apply the block reflector to the
               // residual of the matrix A stored in A(KB+1:M,KB+1:N),
               // since the submatrix is zero and we stop the computation.
               // But, we need to apply the block reflector to the residual
               // right hand sides stored in A(KB+1:M,N+1:N+NRHS), if the
               // residual right hand sides exist.  This occurs
               // when ( NRHS != 0 AND KB <= (M-IOFFSET) ):

               // A(I+1:M,N+1:N+NRHS) := A(I+1:M,N+1:N+NRHS) -
                                // A(I+1:M,1:KB) * F(N+1:N+NRHS,1:KB)**H.

               if ( NRHS.GT.0 .AND. KB.LT.(M-IOFFSET) ) {
                  cgemm('No transpose', 'Conjugate transpose', M-IF, NRHS, KB, -CONE, A( IF+1, 1 ), LDA, F( N+1, 1 ), LDF, CONE, A( IF+1, N+1 ), LDA );
               }

               // There is no need to recompute the 2-norm of the
               // difficult columns, since we stop the factorization.

               // Set TAUs corresponding to the columns that were not
               // factorized to ZERO, i.e. set TAU(KB+1:MINMNFACT) = CZERO,
               // which is equivalent to seting TAU(K:MINMNFACT) = CZERO.

               for (J = K; J <= MINMNFACT; J++) {
                  TAU( J ) = CZERO
               END DO

               // Return from the routine.

               RETURN

            }

            // ============================================================

            // Check if the submatrix A(I:M,K:N) contains Inf,
            // set INFO parameter to the column number, where
            // the first Inf is found plus N, and continue
            // the computation.
            // We need to check the condition only if the
            // column index (same as row index) of the original whole
            // matrix is larger than 1, since the condition for whole
            // original matrix is checked in the main routine.

            if ( INFO.EQ.0 .AND. MAXC2NRMK.GT.HUGEVAL ) {
               INFO = N + K - 1 + KP
            }

            // ============================================================

            // Test for the second and third tolerance stopping criteria.
            // NOTE: There is no need to test for ABSTOL.GE.ZERO, since
            // MAXC2NRMK is non-negative. Similarly, there is no need
            // to test for RELTOL.GE.ZERO, since RELMAXC2NRMK is
            // non-negative.
            // We need to check the condition only if the
            // column index (same as row index) of the original whole
            // matrix is larger than 1, since the condition for whole
            // original matrix is checked in the main routine.

            RELMAXC2NRMK =  MAXC2NRMK / MAXC2NRM

            if ( MAXC2NRMK.LE.ABSTOL .OR. RELMAXC2NRMK.LE.RELTOL ) {

               DONE = .TRUE.

               // Set KB, the number of factorized partial columns
                       // that are non-zero in each step in the block,
                       // i.e. the rank of the factor R.
               // Set IF, the number of processed rows in the block, which
                       // is the same as the number of processed rows in
                       // the original whole matrix A_orig;

                  KB = K - 1
                  IF = I - 1

               // Apply the block reflector to the residual of the
               // matrix A and the residual of the right hand sides B, if
               // the residual matrix and and/or the residual of the right
               // hand sides exist,  i.e. if the submatrix
               // A(I+1:M,KB+1:N+NRHS) exists.  This occurs when
                  // KB < MINMNUPDT = min( M-IOFFSET, N+NRHS ):

               // A(IF+1:M,K+1:N+NRHS) := A(IF+1:M,KB+1:N+NRHS) -
                              // A(IF+1:M,1:KB) * F(KB+1:N+NRHS,1:KB)**H.

               if ( KB.LT.MINMNUPDT ) {
                  cgemm('No transpose', 'Conjugate transpose', M-IF, N+NRHS-KB, KB,-CONE, A( IF+1, 1 ), LDA, F( KB+1, 1 ), LDF, CONE, A( IF+1, KB+1 ), LDA );
               }

               // There is no need to recompute the 2-norm of the
               // difficult columns, since we stop the factorization.

               // Set TAUs corresponding to the columns that were not
               // factorized to ZERO, i.e. set TAU(KB+1:MINMNFACT) = CZERO,
               // which is equivalent to seting TAU(K:MINMNFACT) = CZERO.

               for (J = K; J <= MINMNFACT; J++) {
                  TAU( J ) = CZERO
               END DO

               // Return from the routine.

               RETURN

            }

            // ============================================================

            // End ELSE of IF(I.EQ.1)

         }

         // ===============================================================

         // If the pivot column is not the first column of the
         // subblock A(1:M,K:N):
         // 1) swap the K-th column and the KP-th pivot column
            // in A(1:M,1:N);
         // 2) swap the K-th row and the KP-th row in F(1:N,1:K-1)
         // 3) copy the K-th element into the KP-th element of the partial
            // and exact 2-norm vectors VN1 and VN2. (Swap is not needed
            // for VN1 and VN2 since we use the element with the index
            // larger than K in the next loop step.)
         // 4) Save the pivot interchange with the indices relative to the
            // the original matrix A_orig, not the block A(1:M,1:N).

         if ( KP.NE.K ) {
            cswap(M, A( 1, KP ), 1, A( 1, K ), 1 );
            cswap(K-1, F( KP, 1 ), LDF, F( K, 1 ), LDF );
            VN1( KP ) = VN1( K )
            VN2( KP ) = VN2( K )
            ITEMP = JPIV( KP )
            JPIV( KP ) = JPIV( K )
            JPIV( K ) = ITEMP
         }

         // Apply previous Householder reflectors to column K:
         // A(I:M,K) := A(I:M,K) - A(I:M,1:K-1)*F(K,1:K-1)**H.

         if ( K.GT.1 ) {
            DO J = 1, K - 1
               F( K, J ) = CONJG( F( K, J ) )
            END DO
            cgemv('No transpose', M-I+1, K-1, -CONE, A( I, 1 ), LDA, F( K, 1 ), LDF, CONE, A( I, K ), 1 );
            DO J = 1, K - 1
               F( K, J ) = CONJG( F( K, J ) )
            END DO
         }

         // Generate elementary reflector H(k) using the column A(I:M,K).

         if ( I.LT.M ) {
            clarfg(M-I+1, A( I, K ), A( I+1, K ), 1, TAU( K ) );
         } else {
            TAU( K ) = CZERO
         }

         // Check if TAU(K) contains NaN, set INFO parameter
         // to the column number where NaN is found and return from
         // the routine.
         // NOTE: There is no need to check TAU(K) for Inf,
         // since CLARFG cannot produce TAU(KK) or Householder vector
         // below the diagonal containing Inf. Only BETA on the diagonal,
         // returned by CLARFG can contain Inf, which requires
         // TAU(K) to contain NaN. Therefore, this case of generating Inf
         // by CLARFG is covered by checking TAU(K) for NaN.

         if ( SISNAN( REAL( TAU(K) ) ) ) {
            TAUNAN = REAL( TAU(K) )
         } else if ( SISNAN( AIMAG( TAU(K) ) ) ) {
            TAUNAN = AIMAG( TAU(K) )
         } else {
            TAUNAN = ZERO
         }

         if ( SISNAN( TAUNAN ) ) {

            DONE = .TRUE.

            // Set KB, the number of factorized partial columns
                    // that are non-zero in each step in the block,
                    // i.e. the rank of the factor R.
            // Set IF, the number of processed rows in the block, which
                    // is the same as the number of processed rows in
                    // the original whole matrix A_orig.

            KB = K - 1
            IF = I - 1
            INFO = K

            // Set MAXC2NRMK and  RELMAXC2NRMK to NaN.

            MAXC2NRMK = TAUNAN
            RELMAXC2NRMK = TAUNAN

            // There is no need to apply the block reflector to the
            // residual of the matrix A stored in A(KB+1:M,KB+1:N),
            // since the submatrix contains NaN and we stop
            // the computation.
            // But, we need to apply the block reflector to the residual
            // right hand sides stored in A(KB+1:M,N+1:N+NRHS), if the
            // residual right hand sides exist.  This occurs
            // when ( NRHS != 0 AND KB <= (M-IOFFSET) ):

            // A(I+1:M,N+1:N+NRHS) := A(I+1:M,N+1:N+NRHS) -
                             // A(I+1:M,1:KB) * F(N+1:N+NRHS,1:KB)**H.

            if ( NRHS.GT.0 .AND. KB.LT.(M-IOFFSET) ) {
               cgemm('No transpose', 'Conjugate transpose', M-IF, NRHS, KB, -CONE, A( IF+1, 1 ), LDA, F( N+1, 1 ), LDF, CONE, A( IF+1, N+1 ), LDA );
            }

            // There is no need to recompute the 2-norm of the
            // difficult columns, since we stop the factorization.

            // Array TAU(KF+1:MINMNFACT) is not set and contains
            // undefined elements.

            // Return from the routine.

            RETURN
         }

         // ===============================================================

         AIK = A( I, K )
         A( I, K ) = CONE

         // ===============================================================

         // Compute the current K-th column of F:
           // 1) F(K+1:N,K) := tau(K) * A(I:M,K+1:N)**H * A(I:M,K).

         if ( K.LT.N+NRHS ) {
            cgemv('Conjugate transpose', M-I+1, N+NRHS-K, TAU( K ), A( I, K+1 ), LDA, A( I, K ), 1, CZERO, F( K+1, K ), 1 );
         }

            // 2) Zero out elements above and on the diagonal of the
               // column K in matrix F, i.e elements F(1:K,K).

         for (J = 1; J <= K; J++) {
            F( J, K ) = CZERO
         END DO

          // 3) Incremental updating of the K-th column of F:
         // F(1:N,K) := F(1:N,K) - tau(K) * F(1:N,1:K-1) * A(I:M,1:K-1)**H
                     // * A(I:M,K).

         if ( K.GT.1 ) {
            cgemv('Conjugate Transpose', M-I+1, K-1, -TAU( K ), A( I, 1 ), LDA, A( I, K ), 1, CZERO, AUXV( 1 ), 1 );

            cgemv('No transpose', N+NRHS, K-1, CONE, F( 1, 1 ), LDF, AUXV( 1 ), 1, CONE, F( 1, K ), 1 );
         }

         // ===============================================================

         // Update the current I-th row of A:
         // A(I,K+1:N+NRHS) := A(I,K+1:N+NRHS)
                          // - A(I,1:K)*F(K+1:N+NRHS,1:K)**H.

         if ( K.LT.N+NRHS ) {
            cgemm('No transpose', 'Conjugate transpose', 1, N+NRHS-K, K, -CONE, A( I, 1 ), LDA, F( K+1, 1 ), LDF, CONE, A( I, K+1 ), LDA );
         }

         A( I, K ) = AIK

         // Update the partial column 2-norms for the residual matrix,
         // only if the residual matrix A(I+1:M,K+1:N) exists, i.e.
         // when K < MINMNFACT = min( M-IOFFSET, N ).

         if ( K.LT.MINMNFACT ) {

            DO J = K + 1, N
               if ( VN1( J ).NE.ZERO ) {

                  // NOTE: The following lines follow from the analysis in
                  // Lapack Working Note 176.

                  TEMP = ABS( A( I, J ) ) / VN1( J )
                  TEMP = MAX( ZERO, ( ONE+TEMP )*( ONE-TEMP ) )
                  TEMP2 = TEMP*( VN1( J ) / VN2( J ) )**2
                  if ( TEMP2.LE.TOL3Z ) {

                     // At J-index, we have a difficult column for the
                     // update of the 2-norm. Save the index of the previous
                     // difficult column in IWORK(J-1).
                     // NOTE: ILSTCC > 1, threfore we can use IWORK only
                     // with N-1 elements, where the elements are
                     // shifted by 1 to the left.

                     IWORK( J-1 ) = LSTICC

                     // Set the index of the last difficult column LSTICC.

                     LSTICC = J

                  } else {
                     VN1( J ) = VN1( J )*SQRT( TEMP )
                  }
               }
            END DO

         }

         // End of while loop.

      END DO

      // Now, afler the loop:
         // Set KB, the number of factorized columns in the block;
         // Set IF, the number of processed rows in the block, which
                 // is the same as the number of processed rows in
                 // the original whole matrix A_orig, IF = IOFFSET + KB.

      KB = K
      IF = I

      // Apply the block reflector to the residual of the matrix A
      // and the residual of the right hand sides B, if the residual
      // matrix and and/or the residual of the right hand sides
      // exist,  i.e. if the submatrix A(I+1:M,KB+1:N+NRHS) exists.
      // This occurs when KB < MINMNUPDT = min( M-IOFFSET, N+NRHS ):

      // A(IF+1:M,K+1:N+NRHS) := A(IF+1:M,KB+1:N+NRHS) -
                          // A(IF+1:M,1:KB) * F(KB+1:N+NRHS,1:KB)**H.

      if ( KB.LT.MINMNUPDT ) {
         cgemm('No transpose', 'Conjugate transpose', M-IF, N+NRHS-KB, KB, -CONE, A( IF+1, 1 ), LDA, F( KB+1, 1 ), LDF, CONE, A( IF+1, KB+1 ), LDA );
      }

      // Recompute the 2-norm of the difficult columns.
      // Loop over the index of the difficult columns from the largest
      // to the smallest index.

      DO WHILE( LSTICC.GT.0 )

         // LSTICC is the index of the last difficult column is greater
         // than 1.
         // ITEMP is the index of the previous difficult column.

         ITEMP = IWORK( LSTICC-1 )

         // Compute the 2-norm explicilty for the last difficult column and
         // save it in the partial and exact 2-norm vectors VN1 and VN2.

         // NOTE: The computation of VN1( LSTICC ) relies on the fact that
         // SCNRM2 does not fail on vectors with norm below the value of
         // SQRT(SLAMCH('S'))

         VN1( LSTICC ) = SCNRM2( M-IF, A( IF+1, LSTICC ), 1 )
         VN2( LSTICC ) = VN1( LSTICC )

         // Downdate the index of the last difficult column to
         // the index of the previous difficult column.

         LSTICC = ITEMP

      END DO

      RETURN

      // End of CLAQP3RK

      }