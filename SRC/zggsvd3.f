      SUBROUTINE ZGGSVD3( JOBU, JOBV, JOBQ, M, N, P, K, L, A, LDA, B, LDB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, LWORK, RWORK, IWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBQ, JOBU, JOBV;
      int                INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P, LWORK;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             ALPHA( * ), BETA( * ), RWORK( * );
      COMPLEX*16         A( LDA, * ), B( LDB, * ), Q( LDQ, * ), U( LDU, * ), V( LDV, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               WANTQ, WANTU, WANTV, LQUERY;
      int                I, IBND, ISUB, J, NCYCLE, LWKOPT;
      double             ANORM, BNORM, SMAX, TEMP, TOLA, TOLB, ULP, UNFL;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH, ZLANGE;
      // EXTERNAL LSAME, DLAMCH, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, XERBLA, ZGGSVP3, ZTGSJA
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
      if ( .NOT.( WANTU || LSAME( JOBU, 'N' ) ) ) {
         INFO = -1
      } else if ( .NOT.( WANTV || LSAME( JOBV, 'N' ) ) ) {
         INFO = -2
      } else if ( .NOT.( WANTQ || LSAME( JOBQ, 'N' ) ) ) {
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
      } else if ( LDU.LT.1 || ( WANTU && LDU.LT.M ) ) {
         INFO = -16
      } else if ( LDV.LT.1 || ( WANTV && LDV.LT.P ) ) {
         INFO = -18
      } else if ( LDQ.LT.1 || ( WANTQ && LDQ.LT.N ) ) {
         INFO = -20
      } else if ( LWORK.LT.1 && .NOT.LQUERY ) {
         INFO = -24
      }

      // Compute workspace

      if ( INFO == 0 ) {
         zggsvp3(JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB, TOLA, TOLB, K, L, U, LDU, V, LDV, Q, LDQ, IWORK, RWORK, WORK, WORK, -1, INFO );
         LWKOPT = N + INT( WORK( 1 ) )
         LWKOPT = MAX( 2*N, LWKOPT )
         LWKOPT = MAX( 1, LWKOPT )
         WORK( 1 ) = DCMPLX( LWKOPT )
      }

      if ( INFO != 0 ) {
         xerbla('ZGGSVD3', -INFO );
         RETURN
      }
      if ( LQUERY ) {
         RETURN
      }

      // Compute the Frobenius norm of matrices A and B

      ANORM = ZLANGE( '1', M, N, A, LDA, RWORK )
      BNORM = ZLANGE( '1', P, N, B, LDB, RWORK )

      // Get machine precision and set up threshold for determining
      // the effective numerical rank of the matrices A and B.

      ULP = DLAMCH( 'Precision' )
      UNFL = DLAMCH( 'Safe Minimum' )
      TOLA = MAX( M, N )*MAX( ANORM, UNFL )*ULP
      TOLB = MAX( P, N )*MAX( BNORM, UNFL )*ULP

      zggsvp3(JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB, TOLA, TOLB, K, L, U, LDU, V, LDV, Q, LDQ, IWORK, RWORK, WORK, WORK( N+1 ), LWORK-N, INFO );

      // Compute the GSVD of two upper "triangular" matrices

      ztgsja(JOBU, JOBV, JOBQ, M, P, N, K, L, A, LDA, B, LDB, TOLA, TOLB, ALPHA, BETA, U, LDU, V, LDV, Q, LDQ, WORK, NCYCLE, INFO );

      // Sort the singular values and store the pivot indices in IWORK
      // Copy ALPHA to RWORK, then sort ALPHA in RWORK

      dcopy(N, ALPHA, 1, RWORK, 1 );
      IBND = MIN( L, M-K )
      for (I = 1; I <= IBND; I++) { // 20

         // Scan for largest ALPHA(K+I)

         ISUB = I
         SMAX = RWORK( K+I )
         for (J = I + 1; J <= IBND; J++) { // 10
            TEMP = RWORK( K+J )
            if ( TEMP.GT.SMAX ) {
               ISUB = J
               SMAX = TEMP
            }
         } // 10
         if ( ISUB != I ) {
            RWORK( K+ISUB ) = RWORK( K+I )
            RWORK( K+I ) = SMAX
            IWORK( K+I ) = K + ISUB
         } else {
            IWORK( K+I ) = K + I
         }
      } // 20

      WORK( 1 ) = DCMPLX( LWKOPT )
      RETURN

      // End of ZGGSVD3

      }
