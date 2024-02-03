      void strsyl(TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C, LDC, SCALE, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANA, TRANB;
      int                INFO, ISGN, LDA, LDB, LDC, M, N;
      REAL               SCALE;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), C( LDC, * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               NOTRNA, NOTRNB;
      int                IERR, J, K, K1, K2, KNEXT, L, L1, L2, LNEXT;
      REAL               A11, BIGNUM, DA11, DB, EPS, SCALOC, SGN, SMIN, SMLNUM, SUML, SUMR, XNORM;
      // ..
      // .. Local Arrays ..
      REAL               DUM( 1 ), VEC( 2, 2 ), X( 2, 2 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SDOT, SLAMCH, SLANGE;
      // EXTERNAL LSAME, SDOT, SLAMCH, SLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLALN2, SLASY2, SSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      // Decode and Test input parameters

      NOTRNA = LSAME( TRANA, 'N' );
      NOTRNB = LSAME( TRANB, 'N' );

      INFO = 0;
      if ( !NOTRNA && !LSAME( TRANA, 'T' ) && !LSAME( TRANA, 'C' ) ) {
         INFO = -1;
      } else if ( !NOTRNB && !LSAME( TRANB, 'T' ) && !LSAME( TRANB, 'C' ) ) {
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
         xerbla('STRSYL', -INFO );
         return;
      }

      // Quick return if possible

      SCALE = ONE;
      if (M == 0 || N == 0) return;

      // Set constants to control overflow

      EPS = SLAMCH( 'P' );
      SMLNUM = SLAMCH( 'S' );
      BIGNUM = ONE / SMLNUM;
      SMLNUM = SMLNUM*REAL( M*N ) / EPS;
      BIGNUM = ONE / SMLNUM;

      SMIN = max( SMLNUM, EPS*SLANGE( 'M', M, M, A, LDA, DUM ), EPS*SLANGE( 'M', N, N, B, LDB, DUM ) );

      SGN = ISGN;

      if ( NOTRNA && NOTRNB ) {

         // Solve    A*X + ISGN*X*B = scale*C.

         // The (K,L)th block of X is determined starting from
         // bottom-left corner column by column by

          // A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)

         // Where
                   // M                         L-1
         // R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(J,L)].
                 // I=K+1                       J=1

         // Start column loop (index = L)
         // L1 (L2) : column index of the first (first) row of X(K,L).

         LNEXT = 1;
         for (L = 1; L <= N; L++) { // 70
            if (L < LNEXT) GO TO 70;
            if ( L == N ) {
               L1 = L;
               L2 = L;
            } else {
               if ( B( L+1, L ) != ZERO ) {
                  L1 = L;
                  L2 = L + 1;
                  LNEXT = L + 2;
               } else {
                  L1 = L;
                  L2 = L;
                  LNEXT = L + 1;
               }
            }

            // Start row loop (index = K)
            // K1 (K2): row index of the first (last) row of X(K,L).

            KNEXT = M;
            DO 60 K = M, 1, -1;
               if (K > KNEXT) GO TO 60;
               if ( K == 1 ) {
                  K1 = K;
                  K2 = K;
               } else {
                  if ( A( K, K-1 ) != ZERO ) {
                     K1 = K - 1;
                     K2 = K;
                     KNEXT = K - 2;
                  } else {
                     K1 = K;
                     K2 = K;
                     KNEXT = K - 1;
                  }
               }

               if ( L1 == L2 && K1 == K2 ) {
                  SUML = SDOT( M-K1, A( K1, min( K1+1, M ) ), LDA, C( min( K1+1, M ), L1 ), 1 );
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 );
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR );
                  SCALOC = ONE;

                  A11 = A( K1, K1 ) + SGN*B( L1, L1 );
                  DA11 = ( A11 ).abs();
                  if ( DA11 <= SMIN ) {
                     A11 = SMIN;
                     DA11 = SMIN;
                     INFO = 1;
                  }
                  DB = ( VEC( 1, 1 ) ).abs();
                  if ( DA11 < ONE && DB > ONE ) {
                     if (DB > BIGNUM*DA11) SCALOC = ONE / DB;
                  }
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11;

                  if ( SCALOC != ONE ) {
                     for (J = 1; J <= N; J++) { // 10
                        sscal(M, SCALOC, C( 1, J ), 1 );
                     } // 10
                     SCALE = SCALE*SCALOC;
                  }
                  C( K1, L1 ) = X( 1, 1 );

               } else if ( L1 == L2 && K1 != K2 ) {

                  SUML = SDOT( M-K2, A( K1, min( K2+1, M ) ), LDA, C( min( K2+1, M ), L1 ), 1 );
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 );
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR );

                  SUML = SDOT( M-K2, A( K2, min( K2+1, M ) ), LDA, C( min( K2+1, M ), L1 ), 1 );
                  SUMR = SDOT( L1-1, C( K2, 1 ), LDC, B( 1, L1 ), 1 );
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR );

                  slaln2( false , 2, 1, SMIN, ONE, A( K1, K1 ), LDA, ONE, ONE, VEC, 2, -SGN*B( L1, L1 ), ZERO, X, 2, SCALOC, XNORM, IERR );
                  if (IERR != 0) INFO = 1;

                  if ( SCALOC != ONE ) {
                     for (J = 1; J <= N; J++) { // 20
                        sscal(M, SCALOC, C( 1, J ), 1 );
                     } // 20
                     SCALE = SCALE*SCALOC;
                  }
                  C( K1, L1 ) = X( 1, 1 );
                  C( K2, L1 ) = X( 2, 1 );

               } else if ( L1 != L2 && K1 == K2 ) {

                  SUML = SDOT( M-K1, A( K1, min( K1+1, M ) ), LDA, C( min( K1+1, M ), L1 ), 1 );
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 );
                  VEC( 1, 1 ) = SGN*( C( K1, L1 )-( SUML+SGN*SUMR ) );

                  SUML = SDOT( M-K1, A( K1, min( K1+1, M ) ), LDA, C( min( K1+1, M ), L2 ), 1 );
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L2 ), 1 );
                  VEC( 2, 1 ) = SGN*( C( K1, L2 )-( SUML+SGN*SUMR ) );

                  slaln2( true , 2, 1, SMIN, ONE, B( L1, L1 ), LDB, ONE, ONE, VEC, 2, -SGN*A( K1, K1 ), ZERO, X, 2, SCALOC, XNORM, IERR );
                  if (IERR != 0) INFO = 1;

                  if ( SCALOC != ONE ) {
                     for (J = 1; J <= N; J++) { // 40
                        sscal(M, SCALOC, C( 1, J ), 1 );
                     } // 40
                     SCALE = SCALE*SCALOC;
                  }
                  C( K1, L1 ) = X( 1, 1 );
                  C( K1, L2 ) = X( 2, 1 );

               } else if ( L1 != L2 && K1 != K2 ) {

                  SUML = SDOT( M-K2, A( K1, min( K2+1, M ) ), LDA, C( min( K2+1, M ), L1 ), 1 );
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 );
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR );

                  SUML = SDOT( M-K2, A( K1, min( K2+1, M ) ), LDA, C( min( K2+1, M ), L2 ), 1 );
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L2 ), 1 );
                  VEC( 1, 2 ) = C( K1, L2 ) - ( SUML+SGN*SUMR );

                  SUML = SDOT( M-K2, A( K2, min( K2+1, M ) ), LDA, C( min( K2+1, M ), L1 ), 1 );
                  SUMR = SDOT( L1-1, C( K2, 1 ), LDC, B( 1, L1 ), 1 );
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR );

                  SUML = SDOT( M-K2, A( K2, min( K2+1, M ) ), LDA, C( min( K2+1, M ), L2 ), 1 );
                  SUMR = SDOT( L1-1, C( K2, 1 ), LDC, B( 1, L2 ), 1 );
                  VEC( 2, 2 ) = C( K2, L2 ) - ( SUML+SGN*SUMR );

                  slasy2( false , false , ISGN, 2, 2, A( K1, K1 ), LDA, B( L1, L1 ), LDB, VEC, 2, SCALOC, X, 2, XNORM, IERR );
                  if (IERR != 0) INFO = 1;

                  if ( SCALOC != ONE ) {
                     for (J = 1; J <= N; J++) { // 50
                        sscal(M, SCALOC, C( 1, J ), 1 );
                     } // 50
                     SCALE = SCALE*SCALOC;
                  }
                  C( K1, L1 ) = X( 1, 1 );
                  C( K1, L2 ) = X( 1, 2 );
                  C( K2, L1 ) = X( 2, 1 );
                  C( K2, L2 ) = X( 2, 2 );
               }

            } // 60

         } // 70

      } else if ( !NOTRNA && NOTRNB ) {

         // Solve    A**T *X + ISGN*X*B = scale*C.

         // The (K,L)th block of X is determined starting from
         // upper-left corner column by column by

           // A(K,K)**T*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)

         // Where
                    // K-1                          L-1
           // R(K,L) = SUM [A(I,K)**T*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)]
                    // I=1                          J=1

         // Start column loop (index = L)
         // L1 (L2): column index of the first (last) row of X(K,L)

         LNEXT = 1;
         for (L = 1; L <= N; L++) { // 130
            if (L < LNEXT) GO TO 130;
            if ( L == N ) {
               L1 = L;
               L2 = L;
            } else {
               if ( B( L+1, L ) != ZERO ) {
                  L1 = L;
                  L2 = L + 1;
                  LNEXT = L + 2;
               } else {
                  L1 = L;
                  L2 = L;
                  LNEXT = L + 1;
               }
            }

            // Start row loop (index = K)
            // K1 (K2): row index of the first (last) row of X(K,L)

            KNEXT = 1;
            for (K = 1; K <= M; K++) { // 120
               if (K < KNEXT) GO TO 120;
               if ( K == M ) {
                  K1 = K;
                  K2 = K;
               } else {
                  if ( A( K+1, K ) != ZERO ) {
                     K1 = K;
                     K2 = K + 1;
                     KNEXT = K + 2;
                  } else {
                     K1 = K;
                     K2 = K;
                     KNEXT = K + 1;
                  }
               }

               if ( L1 == L2 && K1 == K2 ) {
                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 );
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 );
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR );
                  SCALOC = ONE;

                  A11 = A( K1, K1 ) + SGN*B( L1, L1 );
                  DA11 = ( A11 ).abs();
                  if ( DA11 <= SMIN ) {
                     A11 = SMIN;
                     DA11 = SMIN;
                     INFO = 1;
                  }
                  DB = ( VEC( 1, 1 ) ).abs();
                  if ( DA11 < ONE && DB > ONE ) {
                     if (DB > BIGNUM*DA11) SCALOC = ONE / DB;
                  }
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11;

                  if ( SCALOC != ONE ) {
                     for (J = 1; J <= N; J++) { // 80
                        sscal(M, SCALOC, C( 1, J ), 1 );
                     } // 80
                     SCALE = SCALE*SCALOC;
                  }
                  C( K1, L1 ) = X( 1, 1 );

               } else if ( L1 == L2 && K1 != K2 ) {

                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 );
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 );
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR );

                  SUML = SDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 );
                  SUMR = SDOT( L1-1, C( K2, 1 ), LDC, B( 1, L1 ), 1 );
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR );

                  slaln2( true , 2, 1, SMIN, ONE, A( K1, K1 ), LDA, ONE, ONE, VEC, 2, -SGN*B( L1, L1 ), ZERO, X, 2, SCALOC, XNORM, IERR );
                  if (IERR != 0) INFO = 1;

                  if ( SCALOC != ONE ) {
                     for (J = 1; J <= N; J++) { // 90
                        sscal(M, SCALOC, C( 1, J ), 1 );
                     } // 90
                     SCALE = SCALE*SCALOC;
                  }
                  C( K1, L1 ) = X( 1, 1 );
                  C( K2, L1 ) = X( 2, 1 );

               } else if ( L1 != L2 && K1 == K2 ) {

                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 );
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 );
                  VEC( 1, 1 ) = SGN*( C( K1, L1 )-( SUML+SGN*SUMR ) );

                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 );
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L2 ), 1 );
                  VEC( 2, 1 ) = SGN*( C( K1, L2 )-( SUML+SGN*SUMR ) );

                  slaln2( true , 2, 1, SMIN, ONE, B( L1, L1 ), LDB, ONE, ONE, VEC, 2, -SGN*A( K1, K1 ), ZERO, X, 2, SCALOC, XNORM, IERR );
                  if (IERR != 0) INFO = 1;

                  if ( SCALOC != ONE ) {
                     for (J = 1; J <= N; J++) { // 100
                        sscal(M, SCALOC, C( 1, J ), 1 );
                     } // 100
                     SCALE = SCALE*SCALOC;
                  }
                  C( K1, L1 ) = X( 1, 1 );
                  C( K1, L2 ) = X( 2, 1 );

               } else if ( L1 != L2 && K1 != K2 ) {

                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 );
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 );
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR );

                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 );
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L2 ), 1 );
                  VEC( 1, 2 ) = C( K1, L2 ) - ( SUML+SGN*SUMR );

                  SUML = SDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 );
                  SUMR = SDOT( L1-1, C( K2, 1 ), LDC, B( 1, L1 ), 1 );
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR );

                  SUML = SDOT( K1-1, A( 1, K2 ), 1, C( 1, L2 ), 1 );
                  SUMR = SDOT( L1-1, C( K2, 1 ), LDC, B( 1, L2 ), 1 );
                  VEC( 2, 2 ) = C( K2, L2 ) - ( SUML+SGN*SUMR );

                  slasy2( true , false , ISGN, 2, 2, A( K1, K1 ), LDA, B( L1, L1 ), LDB, VEC, 2, SCALOC, X, 2, XNORM, IERR );
                  if (IERR != 0) INFO = 1;

                  if ( SCALOC != ONE ) {
                     for (J = 1; J <= N; J++) { // 110
                        sscal(M, SCALOC, C( 1, J ), 1 );
                     } // 110
                     SCALE = SCALE*SCALOC;
                  }
                  C( K1, L1 ) = X( 1, 1 );
                  C( K1, L2 ) = X( 1, 2 );
                  C( K2, L1 ) = X( 2, 1 );
                  C( K2, L2 ) = X( 2, 2 );
               }

            } // 120
         } // 130

      } else if ( !NOTRNA && !NOTRNB ) {

         // Solve    A**T*X + ISGN*X*B**T = scale*C.

         // The (K,L)th block of X is determined starting from
         // top-right corner column by column by

            // A(K,K)**T*X(K,L) + ISGN*X(K,L)*B(L,L)**T = C(K,L) - R(K,L)

         // Where
                      // K-1                            N
             // R(K,L) = SUM [A(I,K)**T*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)**T].
                      // I=1                          J=L+1

         // Start column loop (index = L)
         // L1 (L2): column index of the first (last) row of X(K,L)

         LNEXT = N;
         DO 190 L = N, 1, -1;
            if (L > LNEXT) GO TO 190;
            if ( L == 1 ) {
               L1 = L;
               L2 = L;
            } else {
               if ( B( L, L-1 ) != ZERO ) {
                  L1 = L - 1;
                  L2 = L;
                  LNEXT = L - 2;
               } else {
                  L1 = L;
                  L2 = L;
                  LNEXT = L - 1;
               }
            }

            // Start row loop (index = K)
            // K1 (K2): row index of the first (last) row of X(K,L)

            KNEXT = 1;
            for (K = 1; K <= M; K++) { // 180
               if (K < KNEXT) GO TO 180;
               if ( K == M ) {
                  K1 = K;
                  K2 = K;
               } else {
                  if ( A( K+1, K ) != ZERO ) {
                     K1 = K;
                     K2 = K + 1;
                     KNEXT = K + 2;
                  } else {
                     K1 = K;
                     K2 = K;
                     KNEXT = K + 1;
                  }
               }

               if ( L1 == L2 && K1 == K2 ) {
                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 );
                  SUMR = SDOT( N-L1, C( K1, min( L1+1, N ) ), LDC, B( L1, min( L1+1, N ) ), LDB );
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR );
                  SCALOC = ONE;

                  A11 = A( K1, K1 ) + SGN*B( L1, L1 );
                  DA11 = ( A11 ).abs();
                  if ( DA11 <= SMIN ) {
                     A11 = SMIN;
                     DA11 = SMIN;
                     INFO = 1;
                  }
                  DB = ( VEC( 1, 1 ) ).abs();
                  if ( DA11 < ONE && DB > ONE ) {
                     if (DB > BIGNUM*DA11) SCALOC = ONE / DB;
                  }
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11;

                  if ( SCALOC != ONE ) {
                     for (J = 1; J <= N; J++) { // 140
                        sscal(M, SCALOC, C( 1, J ), 1 );
                     } // 140
                     SCALE = SCALE*SCALOC;
                  }
                  C( K1, L1 ) = X( 1, 1 );

               } else if ( L1 == L2 && K1 != K2 ) {

                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 );
                  SUMR = SDOT( N-L2, C( K1, min( L2+1, N ) ), LDC, B( L1, min( L2+1, N ) ), LDB );
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR );

                  SUML = SDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 );
                  SUMR = SDOT( N-L2, C( K2, min( L2+1, N ) ), LDC, B( L1, min( L2+1, N ) ), LDB );
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR );

                  slaln2( true , 2, 1, SMIN, ONE, A( K1, K1 ), LDA, ONE, ONE, VEC, 2, -SGN*B( L1, L1 ), ZERO, X, 2, SCALOC, XNORM, IERR );
                  if (IERR != 0) INFO = 1;

                  if ( SCALOC != ONE ) {
                     for (J = 1; J <= N; J++) { // 150
                        sscal(M, SCALOC, C( 1, J ), 1 );
                     } // 150
                     SCALE = SCALE*SCALOC;
                  }
                  C( K1, L1 ) = X( 1, 1 );
                  C( K2, L1 ) = X( 2, 1 );

               } else if ( L1 != L2 && K1 == K2 ) {

                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 );
                  SUMR = SDOT( N-L2, C( K1, min( L2+1, N ) ), LDC, B( L1, min( L2+1, N ) ), LDB );
                  VEC( 1, 1 ) = SGN*( C( K1, L1 )-( SUML+SGN*SUMR ) );

                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 );
                  SUMR = SDOT( N-L2, C( K1, min( L2+1, N ) ), LDC, B( L2, min( L2+1, N ) ), LDB );
                  VEC( 2, 1 ) = SGN*( C( K1, L2 )-( SUML+SGN*SUMR ) );

                  slaln2( false , 2, 1, SMIN, ONE, B( L1, L1 ), LDB, ONE, ONE, VEC, 2, -SGN*A( K1, K1 ), ZERO, X, 2, SCALOC, XNORM, IERR );
                  if (IERR != 0) INFO = 1;

                  if ( SCALOC != ONE ) {
                     for (J = 1; J <= N; J++) { // 160
                        sscal(M, SCALOC, C( 1, J ), 1 );
                     } // 160
                     SCALE = SCALE*SCALOC;
                  }
                  C( K1, L1 ) = X( 1, 1 );
                  C( K1, L2 ) = X( 2, 1 );

               } else if ( L1 != L2 && K1 != K2 ) {

                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 );
                  SUMR = SDOT( N-L2, C( K1, min( L2+1, N ) ), LDC, B( L1, min( L2+1, N ) ), LDB );
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR );

                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 );
                  SUMR = SDOT( N-L2, C( K1, min( L2+1, N ) ), LDC, B( L2, min( L2+1, N ) ), LDB );
                  VEC( 1, 2 ) = C( K1, L2 ) - ( SUML+SGN*SUMR );

                  SUML = SDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 );
                  SUMR = SDOT( N-L2, C( K2, min( L2+1, N ) ), LDC, B( L1, min( L2+1, N ) ), LDB );
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR );

                  SUML = SDOT( K1-1, A( 1, K2 ), 1, C( 1, L2 ), 1 );
                  SUMR = SDOT( N-L2, C( K2, min( L2+1, N ) ), LDC, B( L2, min(L2+1, N ) ), LDB );
                  VEC( 2, 2 ) = C( K2, L2 ) - ( SUML+SGN*SUMR );

                  slasy2( true , true , ISGN, 2, 2, A( K1, K1 ), LDA, B( L1, L1 ), LDB, VEC, 2, SCALOC, X, 2, XNORM, IERR );
                  if (IERR != 0) INFO = 1;

                  if ( SCALOC != ONE ) {
                     for (J = 1; J <= N; J++) { // 170
                        sscal(M, SCALOC, C( 1, J ), 1 );
                     } // 170
                     SCALE = SCALE*SCALOC;
                  }
                  C( K1, L1 ) = X( 1, 1 );
                  C( K1, L2 ) = X( 1, 2 );
                  C( K2, L1 ) = X( 2, 1 );
                  C( K2, L2 ) = X( 2, 2 );
               }

            } // 180
         } // 190

      } else if ( NOTRNA && !NOTRNB ) {

         // Solve    A*X + ISGN*X*B**T = scale*C.

         // The (K,L)th block of X is determined starting from
         // bottom-right corner column by column by

             // A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L)**T = C(K,L) - R(K,L)

         // Where
                       // M                          N
             // R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)**T].
                     // I=K+1                      J=L+1

         // Start column loop (index = L)
         // L1 (L2): column index of the first (last) row of X(K,L)

         LNEXT = N;
         DO 250 L = N, 1, -1;
            if (L > LNEXT) GO TO 250;
            if ( L == 1 ) {
               L1 = L;
               L2 = L;
            } else {
               if ( B( L, L-1 ) != ZERO ) {
                  L1 = L - 1;
                  L2 = L;
                  LNEXT = L - 2;
               } else {
                  L1 = L;
                  L2 = L;
                  LNEXT = L - 1;
               }
            }

            // Start row loop (index = K)
            // K1 (K2): row index of the first (last) row of X(K,L)

            KNEXT = M;
            DO 240 K = M, 1, -1;
               if (K > KNEXT) GO TO 240;
               if ( K == 1 ) {
                  K1 = K;
                  K2 = K;
               } else {
                  if ( A( K, K-1 ) != ZERO ) {
                     K1 = K - 1;
                     K2 = K;
                     KNEXT = K - 2;
                  } else {
                     K1 = K;
                     K2 = K;
                     KNEXT = K - 1;
                  }
               }

               if ( L1 == L2 && K1 == K2 ) {
                  SUML = SDOT( M-K1, A( K1, min(K1+1, M ) ), LDA, C( min( K1+1, M ), L1 ), 1 )                   SUMR = SDOT( N-L1, C( K1, min( L1+1, N ) ), LDC, B( L1, min( L1+1, N ) ), LDB );
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR );
                  SCALOC = ONE;

                  A11 = A( K1, K1 ) + SGN*B( L1, L1 );
                  DA11 = ( A11 ).abs();
                  if ( DA11 <= SMIN ) {
                     A11 = SMIN;
                     DA11 = SMIN;
                     INFO = 1;
                  }
                  DB = ( VEC( 1, 1 ) ).abs();
                  if ( DA11 < ONE && DB > ONE ) {
                     if (DB > BIGNUM*DA11) SCALOC = ONE / DB;
                  }
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11;

                  if ( SCALOC != ONE ) {
                     for (J = 1; J <= N; J++) { // 200
                        sscal(M, SCALOC, C( 1, J ), 1 );
                     } // 200
                     SCALE = SCALE*SCALOC;
                  }
                  C( K1, L1 ) = X( 1, 1 );

               } else if ( L1 == L2 && K1 != K2 ) {

                  SUML = SDOT( M-K2, A( K1, min( K2+1, M ) ), LDA, C( min( K2+1, M ), L1 ), 1 )                   SUMR = SDOT( N-L2, C( K1, min( L2+1, N ) ), LDC, B( L1, min( L2+1, N ) ), LDB );
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR );

                  SUML = SDOT( M-K2, A( K2, min( K2+1, M ) ), LDA, C( min( K2+1, M ), L1 ), 1 )                   SUMR = SDOT( N-L2, C( K2, min( L2+1, N ) ), LDC, B( L1, min( L2+1, N ) ), LDB );
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR );

                  slaln2( false , 2, 1, SMIN, ONE, A( K1, K1 ), LDA, ONE, ONE, VEC, 2, -SGN*B( L1, L1 ), ZERO, X, 2, SCALOC, XNORM, IERR );
                  if (IERR != 0) INFO = 1;

                  if ( SCALOC != ONE ) {
                     for (J = 1; J <= N; J++) { // 210
                        sscal(M, SCALOC, C( 1, J ), 1 );
                     } // 210
                     SCALE = SCALE*SCALOC;
                  }
                  C( K1, L1 ) = X( 1, 1 );
                  C( K2, L1 ) = X( 2, 1 );

               } else if ( L1 != L2 && K1 == K2 ) {

                  SUML = SDOT( M-K1, A( K1, min( K1+1, M ) ), LDA, C( min( K1+1, M ), L1 ), 1 )                   SUMR = SDOT( N-L2, C( K1, min( L2+1, N ) ), LDC, B( L1, min( L2+1, N ) ), LDB );
                  VEC( 1, 1 ) = SGN*( C( K1, L1 )-( SUML+SGN*SUMR ) );

                  SUML = SDOT( M-K1, A( K1, min( K1+1, M ) ), LDA, C( min( K1+1, M ), L2 ), 1 )                   SUMR = SDOT( N-L2, C( K1, min( L2+1, N ) ), LDC, B( L2, min( L2+1, N ) ), LDB );
                  VEC( 2, 1 ) = SGN*( C( K1, L2 )-( SUML+SGN*SUMR ) );

                  slaln2( false , 2, 1, SMIN, ONE, B( L1, L1 ), LDB, ONE, ONE, VEC, 2, -SGN*A( K1, K1 ), ZERO, X, 2, SCALOC, XNORM, IERR );
                  if (IERR != 0) INFO = 1;

                  if ( SCALOC != ONE ) {
                     for (J = 1; J <= N; J++) { // 220
                        sscal(M, SCALOC, C( 1, J ), 1 );
                     } // 220
                     SCALE = SCALE*SCALOC;
                  }
                  C( K1, L1 ) = X( 1, 1 );
                  C( K1, L2 ) = X( 2, 1 );

               } else if ( L1 != L2 && K1 != K2 ) {

                  SUML = SDOT( M-K2, A( K1, min( K2+1, M ) ), LDA, C( min( K2+1, M ), L1 ), 1 )                   SUMR = SDOT( N-L2, C( K1, min( L2+1, N ) ), LDC, B( L1, min( L2+1, N ) ), LDB );
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR );

                  SUML = SDOT( M-K2, A( K1, min( K2+1, M ) ), LDA, C( min( K2+1, M ), L2 ), 1 )                   SUMR = SDOT( N-L2, C( K1, min( L2+1, N ) ), LDC, B( L2, min( L2+1, N ) ), LDB );
                  VEC( 1, 2 ) = C( K1, L2 ) - ( SUML+SGN*SUMR );

                  SUML = SDOT( M-K2, A( K2, min( K2+1, M ) ), LDA, C( min( K2+1, M ), L1 ), 1 )                   SUMR = SDOT( N-L2, C( K2, min( L2+1, N ) ), LDC, B( L1, min( L2+1, N ) ), LDB );
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR );

                  SUML = SDOT( M-K2, A( K2, min( K2+1, M ) ), LDA, C( min( K2+1, M ), L2 ), 1 )                   SUMR = SDOT( N-L2, C( K2, min( L2+1, N ) ), LDC, B( L2, min( L2+1, N ) ), LDB );
                  VEC( 2, 2 ) = C( K2, L2 ) - ( SUML+SGN*SUMR );

                  slasy2( false , true , ISGN, 2, 2, A( K1, K1 ), LDA, B( L1, L1 ), LDB, VEC, 2, SCALOC, X, 2, XNORM, IERR );
                  if (IERR != 0) INFO = 1;

                  if ( SCALOC != ONE ) {
                     for (J = 1; J <= N; J++) { // 230
                        sscal(M, SCALOC, C( 1, J ), 1 );
                     } // 230
                     SCALE = SCALE*SCALOC;
                  }
                  C( K1, L1 ) = X( 1, 1 );
                  C( K1, L2 ) = X( 1, 2 );
                  C( K2, L1 ) = X( 2, 1 );
                  C( K2, L2 ) = X( 2, 2 );
               }

            } // 240
         } // 250

      }

      return;
      }
