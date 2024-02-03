      SUBROUTINE SGGSVP3( JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB, TOLA, TOLB, K, L, U, LDU, V, LDV, Q, LDQ, IWORK, TAU, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      IMPLICIT NONE

      // .. Scalar Arguments ..
      String             JOBQ, JOBU, JOBV;
      int                INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P, LWORK;
      REAL               TOLA, TOLB
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL               A( LDA, * ), B( LDB, * ), Q( LDQ, * ), TAU( * ), U( LDU, * ), V( LDV, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               FORWRD, WANTQ, WANTU, WANTV, LQUERY;
      int                I, J, LWKOPT;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      REAL               SROUNDUP_LWORK
      // EXTERNAL SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEQP3, SGEQR2, SGERQ2, SLACPY, SLAPMT, SLASET, SORG2R, SORM2R, SORMR2, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      WANTU = LSAME( JOBU, 'U' )
      WANTV = LSAME( JOBV, 'V' )
      WANTQ = LSAME( JOBQ, 'Q' )
      FORWRD = .TRUE.
      LQUERY = ( LWORK.EQ.-1 )
      LWKOPT = 1

      // Test the input arguments

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
      } else if ( LWORK.LT.1 .AND. .NOT.LQUERY ) {
         INFO = -24
      }

      // Compute workspace

      if ( INFO.EQ.0 ) {
         CALL SGEQP3( P, N, B, LDB, IWORK, TAU, WORK, -1, INFO )
         LWKOPT = INT( WORK ( 1 ) )
         if ( WANTV ) {
            LWKOPT = MAX( LWKOPT, P )
         }
         LWKOPT = MAX( LWKOPT, MIN( N, P ) )
         LWKOPT = MAX( LWKOPT, M )
         if ( WANTQ ) {
            LWKOPT = MAX( LWKOPT, N )
         }
         CALL SGEQP3( M, N, A, LDA, IWORK, TAU, WORK, -1, INFO )
         LWKOPT = MAX( LWKOPT, INT( WORK ( 1 ) ) )
         LWKOPT = MAX( 1, LWKOPT )
         WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'SGGSVP3', -INFO )
         RETURN
      }
      if ( LQUERY ) {
         RETURN
      ENDIF

      // QR with column pivoting of B: B*P = V*( S11 S12 )
                                            // (  0   0  )

      DO 10 I = 1, N
         IWORK( I ) = 0
   10 CONTINUE
      CALL SGEQP3( P, N, B, LDB, IWORK, TAU, WORK, LWORK, INFO )

      // Update A := A*P

      CALL SLAPMT( FORWRD, M, N, A, LDA, IWORK )

      // Determine the effective rank of matrix B.

      L = 0
      DO 20 I = 1, MIN( P, N )
         IF( ABS( B( I, I ) ).GT.TOLB ) L = L + 1
   20 CONTINUE

      if ( WANTV ) {

         // Copy the details of V, and form V.

         CALL SLASET( 'Full', P, P, ZERO, ZERO, V, LDV )
         IF( P.GT.1 ) CALL SLACPY( 'Lower', P-1, N, B( 2, 1 ), LDB, V( 2, 1 ), LDV )
         CALL SORG2R( P, P, MIN( P, N ), V, LDV, TAU, WORK, INFO )
      }

      // Clean up B

      DO 40 J = 1, L - 1
         DO 30 I = J + 1, L
            B( I, J ) = ZERO
   30    CONTINUE
   40 CONTINUE
      IF( P.GT.L ) CALL SLASET( 'Full', P-L, N, ZERO, ZERO, B( L+1, 1 ), LDB )

      if ( WANTQ ) {

         // Set Q = I and Update Q := Q*P

         CALL SLASET( 'Full', N, N, ZERO, ONE, Q, LDQ )
         CALL SLAPMT( FORWRD, N, N, Q, LDQ, IWORK )
      }

      if ( P.GE.L .AND. N.NE.L ) {

         // RQ factorization of (S11 S12): ( S11 S12 ) = ( 0 S12 )*Z

         CALL SGERQ2( L, N, B, LDB, TAU, WORK, INFO )

         // Update A := A*Z**T

         CALL SORMR2( 'Right', 'Transpose', M, N, L, B, LDB, TAU, A, LDA, WORK, INFO )

         if ( WANTQ ) {

            // Update Q := Q*Z**T

            CALL SORMR2( 'Right', 'Transpose', N, N, L, B, LDB, TAU, Q, LDQ, WORK, INFO )
         }

         // Clean up B

         CALL SLASET( 'Full', L, N-L, ZERO, ZERO, B, LDB )
         DO 60 J = N - L + 1, N
            DO 50 I = J - N + L + 1, L
               B( I, J ) = ZERO
   50       CONTINUE
   60    CONTINUE

      }

      // Let              N-L     L
                 // A = ( A11    A12 ) M,

      // then the following does the complete QR decomposition of A11:

               // A11 = U*(  0  T12 )*P1**T
                       // (  0   0  )

      DO 70 I = 1, N - L
         IWORK( I ) = 0
   70 CONTINUE
      CALL SGEQP3( M, N-L, A, LDA, IWORK, TAU, WORK, LWORK, INFO )

      // Determine the effective rank of A11

      K = 0
      DO 80 I = 1, MIN( M, N-L )
         IF( ABS( A( I, I ) ).GT.TOLA ) K = K + 1
   80 CONTINUE

      // Update A12 := U**T*A12, where A12 = A( 1:M, N-L+1:N )

      CALL SORM2R( 'Left', 'Transpose', M, L, MIN( M, N-L ), A, LDA, TAU, A( 1, N-L+1 ), LDA, WORK, INFO )

      if ( WANTU ) {

         // Copy the details of U, and form U

         CALL SLASET( 'Full', M, M, ZERO, ZERO, U, LDU )
         IF( M.GT.1 ) CALL SLACPY( 'Lower', M-1, N-L, A( 2, 1 ), LDA, U( 2, 1 ), LDU )
         CALL SORG2R( M, M, MIN( M, N-L ), U, LDU, TAU, WORK, INFO )
      }

      if ( WANTQ ) {

         // Update Q( 1:N, 1:N-L )  = Q( 1:N, 1:N-L )*P1

         CALL SLAPMT( FORWRD, N, N-L, Q, LDQ, IWORK )
      }

      // Clean up A: set the strictly lower triangular part of
      // A(1:K, 1:K) = 0, and A( K+1:M, 1:N-L ) = 0.

      DO 100 J = 1, K - 1
         DO 90 I = J + 1, K
            A( I, J ) = ZERO
   90    CONTINUE
  100 CONTINUE
      IF( M.GT.K ) CALL SLASET( 'Full', M-K, N-L, ZERO, ZERO, A( K+1, 1 ), LDA )

      if ( N-L.GT.K ) {

         // RQ factorization of ( T11 T12 ) = ( 0 T12 )*Z1

         CALL SGERQ2( K, N-L, A, LDA, TAU, WORK, INFO )

         if ( WANTQ ) {

            // Update Q( 1:N,1:N-L ) = Q( 1:N,1:N-L )*Z1**T

            CALL SORMR2( 'Right', 'Transpose', N, N-L, K, A, LDA, TAU, Q, LDQ, WORK, INFO )
         }

         // Clean up A

         CALL SLASET( 'Full', K, N-L-K, ZERO, ZERO, A, LDA )
         DO 120 J = N - L - K + 1, N - L
            DO 110 I = J - N + L + K + 1, K
               A( I, J ) = ZERO
  110       CONTINUE
  120    CONTINUE

      }

      if ( M.GT.K ) {

         // QR factorization of A( K+1:M,N-L+1:N )

         CALL SGEQR2( M-K, L, A( K+1, N-L+1 ), LDA, TAU, WORK, INFO )

         if ( WANTU ) {

            // Update U(:,K+1:M) := U(:,K+1:M)*U1

            CALL SORM2R( 'Right', 'No transpose', M, M-K, MIN( M-K, L ), A( K+1, N-L+1 ), LDA, TAU, U( 1, K+1 ), LDU, WORK, INFO )
         }

         // Clean up

         DO 140 J = N - L + 1, N
            DO 130 I = J - N + K + L + 1, M
               A( I, J ) = ZERO
  130       CONTINUE
  140    CONTINUE

      }

      WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
      RETURN

      // End of SGGSVP3

      }
