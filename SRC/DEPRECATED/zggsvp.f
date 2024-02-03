      SUBROUTINE ZGGSVP( JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB, TOLA, TOLB, K, L, U, LDU, V, LDV, Q, LDQ, IWORK, RWORK, TAU, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBQ, JOBU, JOBV;
      int                INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P;
      double             TOLA, TOLB;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             RWORK( * );
      COMPLEX*16         A( LDA, * ), B( LDB, * ), Q( LDQ, * ), TAU( * ), U( LDU, * ), V( LDV, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               FORWRD, WANTQ, WANTU, WANTV;
      int                I, J;
      COMPLEX*16         T
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGEQPF, ZGEQR2, ZGERQ2, ZLACPY, ZLAPMT, ZLASET, ZUNG2R, ZUNM2R, ZUNMR2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DIMAG, MAX, MIN
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1( T ) = ABS( DBLE( T ) ) + ABS( DIMAG( T ) )
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      WANTU = LSAME( JOBU, 'U' )
      WANTV = LSAME( JOBV, 'V' )
      WANTQ = LSAME( JOBQ, 'Q' )
      FORWRD = .TRUE.

      INFO = 0
      if ( .NOT.( WANTU .OR. LSAME( JOBU, 'N' ) ) ) {
         INFO = -1
      } else if ( .NOT.( WANTV .OR. LSAME( JOBV, 'N' ) ) ) {
         INFO = -2
      } else if ( .NOT.( WANTQ .OR. LSAME( JOBQ, 'N' ) ) ) {
         INFO = -3
      } else if ( M.LT.0 ) {
         INFO = -4
      } else if ( P.LT.0 ) {
         INFO = -5
      } else if ( N.LT.0 ) {
         INFO = -6
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -8
      } else if ( LDB.LT.MAX( 1, P ) ) {
         INFO = -10
      } else if ( LDU.LT.1 .OR. ( WANTU .AND. LDU.LT.M ) ) {
         INFO = -16
      } else if ( LDV.LT.1 .OR. ( WANTV .AND. LDV.LT.P ) ) {
         INFO = -18
      } else if ( LDQ.LT.1 .OR. ( WANTQ .AND. LDQ.LT.N ) ) {
         INFO = -20
      }
      if ( INFO.NE.0 ) {
         xerbla('ZGGSVP', -INFO );
         RETURN
      }

      // QR with column pivoting of B: B*P = V*( S11 S12 )
                                            // (  0   0  )

      for (I = 1; I <= N; I++) { // 10
         IWORK( I ) = 0
      } // 10
      zgeqpf(P, N, B, LDB, IWORK, TAU, WORK, RWORK, INFO );

      // Update A := A*P

      zlapmt(FORWRD, M, N, A, LDA, IWORK );

      // Determine the effective rank of matrix B.

      L = 0
      DO 20 I = 1, MIN( P, N )
         IF( CABS1( B( I, I ) ).GT.TOLB ) L = L + 1
      } // 20

      if ( WANTV ) {

         // Copy the details of V, and form V.

         zlaset('Full', P, P, CZERO, CZERO, V, LDV );
         IF( P.GT.1 ) CALL ZLACPY( 'Lower', P-1, N, B( 2, 1 ), LDB, V( 2, 1 ), LDV )
         zung2r(P, P, MIN( P, N ), V, LDV, TAU, WORK, INFO );
      }

      // Clean up B

      DO 40 J = 1, L - 1
         DO 30 I = J + 1, L
            B( I, J ) = CZERO
         } // 30
      } // 40
      IF( P.GT.L ) CALL ZLASET( 'Full', P-L, N, CZERO, CZERO, B( L+1, 1 ), LDB )

      if ( WANTQ ) {

         // Set Q = I and Update Q := Q*P

         zlaset('Full', N, N, CZERO, CONE, Q, LDQ );
         zlapmt(FORWRD, N, N, Q, LDQ, IWORK );
      }

      if ( P.GE.L .AND. N.NE.L ) {

         // RQ factorization of ( S11 S12 ) = ( 0 S12 )*Z

         zgerq2(L, N, B, LDB, TAU, WORK, INFO );

         // Update A := A*Z**H

         zunmr2('Right', 'Conjugate transpose', M, N, L, B, LDB, TAU, A, LDA, WORK, INFO );
         if ( WANTQ ) {

            // Update Q := Q*Z**H

            zunmr2('Right', 'Conjugate transpose', N, N, L, B, LDB, TAU, Q, LDQ, WORK, INFO );
         }

         // Clean up B

         zlaset('Full', L, N-L, CZERO, CZERO, B, LDB );
         DO 60 J = N - L + 1, N
            DO 50 I = J - N + L + 1, L
               B( I, J ) = CZERO
            } // 50
         } // 60

      }

      // Let              N-L     L
                 // A = ( A11    A12 ) M,

      // then the following does the complete QR decomposition of A11:

               // A11 = U*(  0  T12 )*P1**H
                       // (  0   0  )

      DO 70 I = 1, N - L
         IWORK( I ) = 0
      } // 70
      zgeqpf(M, N-L, A, LDA, IWORK, TAU, WORK, RWORK, INFO );

      // Determine the effective rank of A11

      K = 0
      DO 80 I = 1, MIN( M, N-L )
         IF( CABS1( A( I, I ) ).GT.TOLA ) K = K + 1
      } // 80

      // Update A12 := U**H*A12, where A12 = A( 1:M, N-L+1:N )

      zunm2r('Left', 'Conjugate transpose', M, L, MIN( M, N-L ), A, LDA, TAU, A( 1, N-L+1 ), LDA, WORK, INFO );

      if ( WANTU ) {

         // Copy the details of U, and form U

         zlaset('Full', M, M, CZERO, CZERO, U, LDU );
         IF( M.GT.1 ) CALL ZLACPY( 'Lower', M-1, N-L, A( 2, 1 ), LDA, U( 2, 1 ), LDU )
         zung2r(M, M, MIN( M, N-L ), U, LDU, TAU, WORK, INFO );
      }

      if ( WANTQ ) {

         // Update Q( 1:N, 1:N-L )  = Q( 1:N, 1:N-L )*P1

         zlapmt(FORWRD, N, N-L, Q, LDQ, IWORK );
      }

      // Clean up A: set the strictly lower triangular part of
      // A(1:K, 1:K) = 0, and A( K+1:M, 1:N-L ) = 0.

      DO 100 J = 1, K - 1
         DO 90 I = J + 1, K
            A( I, J ) = CZERO
         } // 90
      } // 100
      IF( M.GT.K ) CALL ZLASET( 'Full', M-K, N-L, CZERO, CZERO, A( K+1, 1 ), LDA )

      if ( N-L.GT.K ) {

         // RQ factorization of ( T11 T12 ) = ( 0 T12 )*Z1

         zgerq2(K, N-L, A, LDA, TAU, WORK, INFO );

         if ( WANTQ ) {

            // Update Q( 1:N,1:N-L ) = Q( 1:N,1:N-L )*Z1**H

            zunmr2('Right', 'Conjugate transpose', N, N-L, K, A, LDA, TAU, Q, LDQ, WORK, INFO );
         }

         // Clean up A

         zlaset('Full', K, N-L-K, CZERO, CZERO, A, LDA );
         DO 120 J = N - L - K + 1, N - L
            DO 110 I = J - N + L + K + 1, K
               A( I, J ) = CZERO
            } // 110
         } // 120

      }

      if ( M.GT.K ) {

         // QR factorization of A( K+1:M,N-L+1:N )

         zgeqr2(M-K, L, A( K+1, N-L+1 ), LDA, TAU, WORK, INFO );

         if ( WANTU ) {

            // Update U(:,K+1:M) := U(:,K+1:M)*U1

            zunm2r('Right', 'No transpose', M, M-K, MIN( M-K, L ), A( K+1, N-L+1 ), LDA, TAU, U( 1, K+1 ), LDU, WORK, INFO );
         }

         // Clean up

         DO 140 J = N - L + 1, N
            DO 130 I = J - N + K + L + 1, M
               A( I, J ) = CZERO
            } // 130
         } // 140

      }

      RETURN

      // End of ZGGSVP

      }
