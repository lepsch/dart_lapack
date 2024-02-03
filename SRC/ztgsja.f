      SUBROUTINE ZTGSJA( JOBU, JOBV, JOBQ, M, P, N, K, L, A, LDA, B, LDB, TOLA, TOLB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, NCYCLE, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBQ, JOBU, JOBV;
      int                INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, NCYCLE, P;
      double             TOLA, TOLB;
      // ..
      // .. Array Arguments ..
      double             ALPHA( * ), BETA( * );
      COMPLEX*16         A( LDA, * ), B( LDB, * ), Q( LDQ, * ), U( LDU, * ), V( LDV, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      int                MAXIT;
      const              MAXIT = 40 ;
      double             ZERO, ONE, HUGENUM;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..

      bool               INITQ, INITU, INITV, UPPER, WANTQ, WANTU, WANTV;
      int                I, J, KCYCLE;
      double             A1, A3, B1, B3, CSQ, CSU, CSV, ERROR, GAMMA, RWK, SSMIN;
      COMPLEX*16         A2, B2, SNQ, SNU, SNV
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLARTG, XERBLA, ZCOPY, ZDSCAL, ZLAGS2, ZLAPLL, ZLASET, ZROT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCONJG, MAX, MIN, HUGE
      const              HUGENUM = HUGE(ZERO) ;
      // ..
      // .. Executable Statements ..

      // Decode and test the input parameters

      INITU = LSAME( JOBU, 'I' )
      WANTU = INITU .OR. LSAME( JOBU, 'U' )

      INITV = LSAME( JOBV, 'I' )
      WANTV = INITV .OR. LSAME( JOBV, 'V' )

      INITQ = LSAME( JOBQ, 'I' )
      WANTQ = INITQ .OR. LSAME( JOBQ, 'Q' )

      INFO = 0
      if ( .NOT.( INITU .OR. WANTU .OR. LSAME( JOBU, 'N' ) ) ) {
         INFO = -1
      } else if ( .NOT.( INITV .OR. WANTV .OR. LSAME( JOBV, 'N' ) ) ) {
         INFO = -2
      } else if ( .NOT.( INITQ .OR. WANTQ .OR. LSAME( JOBQ, 'N' ) ) ) {
         INFO = -3
      } else if ( M.LT.0 ) {
         INFO = -4
      } else if ( P.LT.0 ) {
         INFO = -5
      } else if ( N.LT.0 ) {
         INFO = -6
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -10
      } else if ( LDB.LT.MAX( 1, P ) ) {
         INFO = -12
      } else if ( LDU.LT.1 .OR. ( WANTU .AND. LDU.LT.M ) ) {
         INFO = -18
      } else if ( LDV.LT.1 .OR. ( WANTV .AND. LDV.LT.P ) ) {
         INFO = -20
      } else if ( LDQ.LT.1 .OR. ( WANTQ .AND. LDQ.LT.N ) ) {
         INFO = -22
      }
      if ( INFO.NE.0 ) {
         xerbla('ZTGSJA', -INFO );
         RETURN
      }

      // Initialize U, V and Q, if necessary

      IF( INITU ) CALL ZLASET( 'Full', M, M, CZERO, CONE, U, LDU )       IF( INITV ) CALL ZLASET( 'Full', P, P, CZERO, CONE, V, LDV )       IF( INITQ ) CALL ZLASET( 'Full', N, N, CZERO, CONE, Q, LDQ )

      // Loop until convergence

      UPPER = .FALSE.
      for (KCYCLE = 1; KCYCLE <= MAXIT; KCYCLE++) { // 40

         UPPER = .NOT.UPPER

         DO 20 I = 1, L - 1
            DO 10 J = I + 1, L

               A1 = ZERO
               A2 = CZERO
               A3 = ZERO
               IF( K+I.LE.M ) A1 = DBLE( A( K+I, N-L+I ) )                IF( K+J.LE.M ) A3 = DBLE( A( K+J, N-L+J ) )

               B1 = DBLE( B( I, N-L+I ) )
               B3 = DBLE( B( J, N-L+J ) )

               if ( UPPER ) {
                  IF( K+I.LE.M ) A2 = A( K+I, N-L+J )
                  B2 = B( I, N-L+J )
               } else {
                  IF( K+J.LE.M ) A2 = A( K+J, N-L+I )
                  B2 = B( J, N-L+I )
               }

               zlags2(UPPER, A1, A2, A3, B1, B2, B3, CSU, SNU, CSV, SNV, CSQ, SNQ );

               // Update (K+I)-th and (K+J)-th rows of matrix A: U**H *A

               IF( K+J.LE.M ) CALL ZROT( L, A( K+J, N-L+1 ), LDA, A( K+I, N-L+1 ), LDA, CSU, DCONJG( SNU ) )

               // Update I-th and J-th rows of matrix B: V**H *B

               zrot(L, B( J, N-L+1 ), LDB, B( I, N-L+1 ), LDB, CSV, DCONJG( SNV ) );

               // Update (N-L+I)-th and (N-L+J)-th columns of matrices
               // A and B: A*Q and B*Q

               zrot(MIN( K+L, M ), A( 1, N-L+J ), 1, A( 1, N-L+I ), 1, CSQ, SNQ );

               zrot(L, B( 1, N-L+J ), 1, B( 1, N-L+I ), 1, CSQ, SNQ );

               if ( UPPER ) {
                  IF( K+I.LE.M ) A( K+I, N-L+J ) = CZERO
                  B( I, N-L+J ) = CZERO
               } else {
                  IF( K+J.LE.M ) A( K+J, N-L+I ) = CZERO
                  B( J, N-L+I ) = CZERO
               }

               // Ensure that the diagonal elements of A and B are real.

               IF( K+I.LE.M ) A( K+I, N-L+I ) = DBLE( A( K+I, N-L+I ) )                IF( K+J.LE.M ) A( K+J, N-L+J ) = DBLE( A( K+J, N-L+J ) )
               B( I, N-L+I ) = DBLE( B( I, N-L+I ) )
               B( J, N-L+J ) = DBLE( B( J, N-L+J ) )

               // Update unitary matrices U, V, Q, if desired.

               IF( WANTU .AND. K+J.LE.M ) CALL ZROT( M, U( 1, K+J ), 1, U( 1, K+I ), 1, CSU, SNU )

               IF( WANTV ) CALL ZROT( P, V( 1, J ), 1, V( 1, I ), 1, CSV, SNV )

               IF( WANTQ ) CALL ZROT( N, Q( 1, N-L+J ), 1, Q( 1, N-L+I ), 1, CSQ, SNQ )

            } // 10
         } // 20

         if ( .NOT.UPPER ) {

            // The matrices A13 and B13 were lower triangular at the start
            // of the cycle, and are now upper triangular.

            // Convergence test: test the parallelism of the corresponding
            // rows of A and B.

            ERROR = ZERO
            DO 30 I = 1, MIN( L, M-K )
               zcopy(L-I+1, A( K+I, N-L+I ), LDA, WORK, 1 );
               zcopy(L-I+1, B( I, N-L+I ), LDB, WORK( L+1 ), 1 );
               zlapll(L-I+1, WORK, 1, WORK( L+1 ), 1, SSMIN );
               ERROR = MAX( ERROR, SSMIN )
            } // 30

            IF( ABS( ERROR ).LE.MIN( TOLA, TOLB ) ) GO TO 50
         }

         // End of cycle loop

      } // 40

      // The algorithm has not converged after MAXIT cycles.

      INFO = 1
      GO TO 100

      } // 50

      // If ERROR <= MIN(TOLA,TOLB), then the algorithm has converged.
      // Compute the generalized singular value pairs (ALPHA, BETA), and
      // set the triangular matrix R to array A.

      for (I = 1; I <= K; I++) { // 60
         ALPHA( I ) = ONE
         BETA( I ) = ZERO
      } // 60

      DO 70 I = 1, MIN( L, M-K )

         A1 = DBLE( A( K+I, N-L+I ) )
         B1 = DBLE( B( I, N-L+I ) )
         GAMMA = B1 / A1

         if ( (GAMMA.LE.HUGENUM).AND.(GAMMA.GE.-HUGENUM) ) {

            if ( GAMMA.LT.ZERO ) {
               zdscal(L-I+1, -ONE, B( I, N-L+I ), LDB );
               IF( WANTV ) CALL ZDSCAL( P, -ONE, V( 1, I ), 1 )
            }

            dlartg(ABS( GAMMA ), ONE, BETA( K+I ), ALPHA( K+I ), RWK );

            if ( ALPHA( K+I ).GE.BETA( K+I ) ) {
               zdscal(L-I+1, ONE / ALPHA( K+I ), A( K+I, N-L+I ), LDA );
            } else {
               zdscal(L-I+1, ONE / BETA( K+I ), B( I, N-L+I ), LDB )                CALL ZCOPY( L-I+1, B( I, N-L+I ), LDB, A( K+I, N-L+I ), LDA );
            }

         } else {

            ALPHA( K+I ) = ZERO
            BETA( K+I ) = ONE
            zcopy(L-I+1, B( I, N-L+I ), LDB, A( K+I, N-L+I ), LDA );
         }
      } // 70

      // Post-assignment

      DO 80 I = M + 1, K + L
         ALPHA( I ) = ZERO
         BETA( I ) = ONE
      } // 80

      if ( K+L.LT.N ) {
         DO 90 I = K + L + 1, N
            ALPHA( I ) = ZERO
            BETA( I ) = ZERO
         } // 90
      }

      } // 100
      NCYCLE = KCYCLE

      RETURN

      // End of ZTGSJA

      }
