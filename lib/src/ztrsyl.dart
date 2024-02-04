      void ztrsyl(TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C, LDC, SCALE, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANA, TRANB;
      int                INFO, ISGN, LDA, LDB, LDC, M, N;
      double             SCALE;
      // ..
      // .. Array Arguments ..
      Complex         A( LDA, * ), B( LDB, * ), C( LDC, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE;
      const              ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               NOTRNA, NOTRNB;
      int                J, K, L;
      double             BIGNUM, DA11, DB, EPS, SCALOC, SGN, SMIN, SMLNUM;
      Complex         A11, SUML, SUMR, VEC, X11;
      // ..
      // .. Local Arrays ..
      double             DUM( 1 );
      // ..
      // .. External Functions ..
      //- bool               LSAME;
      //- double             DLAMCH, ZLANGE;
      //- Complex         ZDOTC, ZDOTU, ZLADIV;
      // EXTERNAL LSAME, DLAMCH, ZLANGE, ZDOTC, ZDOTU, ZLADIV
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZDSCAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCMPLX, DCONJG, DIMAG, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Decode and Test input parameters

      NOTRNA = LSAME( TRANA, 'N' );
      NOTRNB = LSAME( TRANB, 'N' );

      INFO = 0;
      if ( !NOTRNA && !LSAME( TRANA, 'C' ) ) {
         INFO = -1;
      } else if ( !NOTRNB && !LSAME( TRANB, 'C' ) ) {
         INFO = -2;
      } else if ( ISGN != 1 && ISGN != -1 ) {
         INFO = -3;
      } else if ( M < 0 ) {
         INFO = -4;
      } else if ( N < 0 ) {
         INFO = -5;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -7;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -9;
      } else if ( LDC < max( 1, M ) ) {
         INFO = -11;
      }
      if ( INFO != 0 ) {
         xerbla('ZTRSYL', -INFO );
         return;
      }

      // Quick return if possible

      SCALE = ONE;
      if (M == 0 || N == 0) return;

      // Set constants to control overflow

      EPS = DLAMCH( 'P' );
      SMLNUM = DLAMCH( 'S' );
      BIGNUM = ONE / SMLNUM;
      SMLNUM = SMLNUM*(M*N).toDouble() / EPS;
      BIGNUM = ONE / SMLNUM;
      SMIN = max( SMLNUM, EPS*ZLANGE( 'M', M, M, A, LDA, DUM ), EPS*ZLANGE( 'M', N, N, B, LDB, DUM ) );
      SGN = ISGN;

      if ( NOTRNA && NOTRNB ) {

         // Solve    A*X + ISGN*X*B = scale*C.

         // The (K,L)th block of X is determined starting from
         // bottom-left corner column by column by

             // A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)

         // Where
                     // M                        L-1
           // R(K,L) = SUM [A(K,I)*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)].
                   // I=K+1                      J=1

         for (L = 1; L <= N; L++) { // 30
            for (K = M; K >= 1; K--) { // 20

               SUML = ZDOTU( M-K, A( K, min( K+1, M ) ), LDA, C( min( K+1, M ), L ), 1 );
               SUMR = ZDOTU( L-1, C( K, 1 ), LDC, B( 1, L ), 1 );
               VEC = C( K, L ) - ( SUML+SGN*SUMR );

               SCALOC = ONE;
               A11 = A( K, K ) + SGN*B( L, L );
               DA11 = ( A11.toDouble() ).abs() + ( DIMAG( A11 ) ).abs();
               if ( DA11 <= SMIN ) {
                  A11 = SMIN;
                  DA11 = SMIN;
                  INFO = 1;
               }
               DB = ( VEC.toDouble() ).abs() + ( DIMAG( VEC ) ).abs();
               if ( DA11 < ONE && DB > ONE ) {
                  if (DB > BIGNUM*DA11) SCALOC = ONE / DB;
               }
               X11 = ZLADIV( VEC*DCMPLX( SCALOC ), A11 );

               if ( SCALOC != ONE ) {
                  for (J = 1; J <= N; J++) { // 10
                     zdscal(M, SCALOC, C( 1, J ), 1 );
                  } // 10
                  SCALE = SCALE*SCALOC;
               }
               C[K, L] = X11;

            } // 20
         } // 30

      } else if ( !NOTRNA && NOTRNB ) {

         // Solve    A**H *X + ISGN*X*B = scale*C.

         // The (K,L)th block of X is determined starting from
         // upper-left corner column by column by

             // A**H(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)

         // Where
                    // K-1                           L-1
           // R(K,L) = SUM [A**H(I,K)*X(I,L)] + ISGN*SUM [X(K,J)*B(J,L)]
                    // I=1                           J=1

         for (L = 1; L <= N; L++) { // 60
            for (K = 1; K <= M; K++) { // 50

               SUML = ZDOTC( K-1, A( 1, K ), 1, C( 1, L ), 1 );
               SUMR = ZDOTU( L-1, C( K, 1 ), LDC, B( 1, L ), 1 );
               VEC = C( K, L ) - ( SUML+SGN*SUMR );

               SCALOC = ONE;
               A11 = DCONJG( A( K, K ) ) + SGN*B( L, L );
               DA11 = ( A11.toDouble() ).abs() + ( DIMAG( A11 ) ).abs();
               if ( DA11 <= SMIN ) {
                  A11 = SMIN;
                  DA11 = SMIN;
                  INFO = 1;
               }
               DB = ( VEC.toDouble() ).abs() + ( DIMAG( VEC ) ).abs();
               if ( DA11 < ONE && DB > ONE ) {
                  if (DB > BIGNUM*DA11) SCALOC = ONE / DB;
               }

               X11 = ZLADIV( VEC*DCMPLX( SCALOC ), A11 );

               if ( SCALOC != ONE ) {
                  for (J = 1; J <= N; J++) { // 40
                     zdscal(M, SCALOC, C( 1, J ), 1 );
                  } // 40
                  SCALE = SCALE*SCALOC;
               }
               C[K, L] = X11;

            } // 50
         } // 60

      } else if ( !NOTRNA && !NOTRNB ) {

         // Solve    A**H*X + ISGN*X*B**H = C.

         // The (K,L)th block of X is determined starting from
         // upper-right corner column by column by

             // A**H(K,K)*X(K,L) + ISGN*X(K,L)*B**H(L,L) = C(K,L) - R(K,L)

         // Where
                     // K-1
            // R(K,L) = SUM [A**H(I,K)*X(I,L)] +
                     // I=1
                            // N
                      // ISGN*SUM [X(K,J)*B**H(L,J)].
                           // J=L+1

         for (L = N; L >= 1; L--) { // 90
            for (K = 1; K <= M; K++) { // 80

               SUML = ZDOTC( K-1, A( 1, K ), 1, C( 1, L ), 1 );
               SUMR = ZDOTC( N-L, C( K, min( L+1, N ) ), LDC, B( L, min( L+1, N ) ), LDB );
               VEC = C( K, L ) - ( SUML+SGN*DCONJG( SUMR ) );

               SCALOC = ONE;
               A11 = DCONJG( A( K, K )+SGN*B( L, L ) );
               DA11 = ( A11.toDouble() ).abs() + ( DIMAG( A11 ) ).abs();
               if ( DA11 <= SMIN ) {
                  A11 = SMIN;
                  DA11 = SMIN;
                  INFO = 1;
               }
               DB = ( VEC.toDouble() ).abs() + ( DIMAG( VEC ) ).abs();
               if ( DA11 < ONE && DB > ONE ) {
                  if (DB > BIGNUM*DA11) SCALOC = ONE / DB;
               }

               X11 = ZLADIV( VEC*DCMPLX( SCALOC ), A11 );

               if ( SCALOC != ONE ) {
                  for (J = 1; J <= N; J++) { // 70
                     zdscal(M, SCALOC, C( 1, J ), 1 );
                  } // 70
                  SCALE = SCALE*SCALOC;
               }
               C[K, L] = X11;

            } // 80
         } // 90

      } else if ( NOTRNA && !NOTRNB ) {

         // Solve    A*X + ISGN*X*B**H = C.

         // The (K,L)th block of X is determined starting from
         // bottom-left corner column by column by

            // A(K,K)*X(K,L) + ISGN*X(K,L)*B**H(L,L) = C(K,L) - R(K,L)

         // Where
                     // M                          N
           // R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B**H(L,J)]
                   // I=K+1                      J=L+1

         for (L = N; L >= 1; L--) { // 120
            for (K = M; K >= 1; K--) { // 110

               SUML = ZDOTU( M-K, A( K, min( K+1, M ) ), LDA, C( min( K+1, M ), L ), 1 )                SUMR = ZDOTC( N-L, C( K, min( L+1, N ) ), LDC, B( L, min( L+1, N ) ), LDB );
               VEC = C( K, L ) - ( SUML+SGN*DCONJG( SUMR ) );

               SCALOC = ONE;
               A11 = A( K, K ) + SGN*DCONJG( B( L, L ) );
               DA11 = ( A11.toDouble() ).abs() + ( DIMAG( A11 ) ).abs();
               if ( DA11 <= SMIN ) {
                  A11 = SMIN;
                  DA11 = SMIN;
                  INFO = 1;
               }
               DB = ( VEC.toDouble() ).abs() + ( DIMAG( VEC ) ).abs();
               if ( DA11 < ONE && DB > ONE ) {
                  if (DB > BIGNUM*DA11) SCALOC = ONE / DB;
               }

               X11 = ZLADIV( VEC*DCMPLX( SCALOC ), A11 );

               if ( SCALOC != ONE ) {
                  for (J = 1; J <= N; J++) { // 100
                     zdscal(M, SCALOC, C( 1, J ), 1 );
                  } // 100
                  SCALE = SCALE*SCALOC;
               }
               C[K, L] = X11;

            } // 110
         } // 120

      }

      return;
      }
