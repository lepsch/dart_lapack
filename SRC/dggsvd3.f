      SUBROUTINE DGGSVD3( JOBU, JOBV, JOBQ, M, N, P, K, L, A, LDA, B, LDB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, LWORK, IWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBQ, JOBU, JOBV;
      int                INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P, LWORK;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             A( LDA, * ), ALPHA( * ), B( LDB, * ), BETA( * ), Q( LDQ, * ), U( LDU, * ), V( LDV, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               WANTQ, WANTU, WANTV, LQUERY;
      int                I, IBND, ISUB, J, NCYCLE, LWKOPT;
      double             ANORM, BNORM, SMAX, TEMP, TOLA, TOLB, ULP, UNFL;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, DLANGE;
      // EXTERNAL LSAME, DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DGGSVP3, DTGSJA, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Decode and test the input parameters

      WANTU = LSAME( JOBU, 'U' )
      WANTV = LSAME( JOBV, 'V' )
      WANTQ = LSAME( JOBQ, 'Q' )
      LQUERY = ( LWORK == -1 )
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
      } else if ( N.LT.0 ) {
         INFO = -5
      } else if ( P.LT.0 ) {
         INFO = -6
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -10
      } else if ( LDB.LT.MAX( 1, P ) ) {
         INFO = -12
      } else if ( LDU.LT.1 .OR. ( WANTU && LDU.LT.M ) ) {
         INFO = -16
      } else if ( LDV.LT.1 .OR. ( WANTV && LDV.LT.P ) ) {
         INFO = -18
      } else if ( LDQ.LT.1 .OR. ( WANTQ && LDQ.LT.N ) ) {
         INFO = -20
      } else if ( LWORK.LT.1 && .NOT.LQUERY ) {
         INFO = -24
      }

      // Compute workspace

      if ( INFO == 0 ) {
         dggsvp3(JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB, TOLA, TOLB, K, L, U, LDU, V, LDV, Q, LDQ, IWORK, WORK, WORK, -1, INFO );
         LWKOPT = N + INT( WORK( 1 ) )
         LWKOPT = MAX( 2*N, LWKOPT )
         LWKOPT = MAX( 1, LWKOPT )
         WORK( 1 ) = DBLE( LWKOPT )
      }

      if ( INFO != 0 ) {
         xerbla('DGGSVD3', -INFO );
         RETURN
      }
      if ( LQUERY ) {
         RETURN
      }

      // Compute the Frobenius norm of matrices A and B

      ANORM = DLANGE( '1', M, N, A, LDA, WORK )
      BNORM = DLANGE( '1', P, N, B, LDB, WORK )

      // Get machine precision and set up threshold for determining
      // the effective numerical rank of the matrices A and B.

      ULP = DLAMCH( 'Precision' )
      UNFL = DLAMCH( 'Safe Minimum' )
      TOLA = MAX( M, N )*MAX( ANORM, UNFL )*ULP
      TOLB = MAX( P, N )*MAX( BNORM, UNFL )*ULP

      // Preprocessing

      dggsvp3(JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB, TOLA, TOLB, K, L, U, LDU, V, LDV, Q, LDQ, IWORK, WORK, WORK( N+1 ), LWORK-N, INFO );

      // Compute the GSVD of two upper "triangular" matrices

      dtgsja(JOBU, JOBV, JOBQ, M, P, N, K, L, A, LDA, B, LDB, TOLA, TOLB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, NCYCLE, INFO );

      // Sort the singular values and store the pivot indices in IWORK
      // Copy ALPHA to WORK, then sort ALPHA in WORK

      dcopy(N, ALPHA, 1, WORK, 1 );
      IBND = MIN( L, M-K )
      for (I = 1; I <= IBND; I++) { // 20

         // Scan for largest ALPHA(K+I)

         ISUB = I
         SMAX = WORK( K+I )
         for (J = I + 1; J <= IBND; J++) { // 10
            TEMP = WORK( K+J )
            if ( TEMP.GT.SMAX ) {
               ISUB = J
               SMAX = TEMP
            }
         } // 10
         if ( ISUB != I ) {
            WORK( K+ISUB ) = WORK( K+I )
            WORK( K+I ) = SMAX
            IWORK( K+I ) = K + ISUB
         } else {
            IWORK( K+I ) = K + I
         }
      } // 20

      WORK( 1 ) = DBLE( LWKOPT )
      RETURN

      // End of DGGSVD3

      }
