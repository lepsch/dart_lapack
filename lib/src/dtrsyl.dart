import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dtrsyl(TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C, LDC, SCALE, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANA, TRANB;
      int                INFO, ISGN, LDA, LDB, LDC, M, N;
      double             SCALE;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), B( LDB, * ), C( LDC, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               NOTRNA, NOTRNB;
      int                IERR, J, K, K1, K2, KNEXT, L, L1, L2, LNEXT;
      double             A11, BIGNUM, DA11, DB, EPS, SCALOC, SGN, SMIN, SMLNUM, SUML, SUMR, XNORM;
      // ..
      // .. Local Arrays ..
      double             DUM( 1 ), VEC( 2, 2 ), X( 2, 2 );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DDOT, DLAMCH, DLANGE;
      // EXTERNAL lsame, DDOT, DLAMCH, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLALN2, DLASY2, DSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Decode and Test input parameters

      NOTRNA = lsame( TRANA, 'N' );
      NOTRNB = lsame( TRANB, 'N' );

      INFO = 0;
      if ( !NOTRNA && !lsame( TRANA, 'T' ) && !lsame( TRANA, 'C' ) ) {
         INFO = -1;
      } else if ( !NOTRNB && !lsame( TRANB, 'T' ) && !lsame( TRANB, 'C' ) ) {
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
         xerbla('DTRSYL', -INFO );
         return;
      }

      // Quick return if possible

      SCALE = ONE;
      if (M == 0 || N == 0) return;

      // Set constants to control overflow

      EPS = dlamch( 'P' );
      SMLNUM = dlamch( 'S' );
      BIGNUM = ONE / SMLNUM;
      SMLNUM = SMLNUM*(M*N).toDouble() / EPS;
      BIGNUM = ONE / SMLNUM;

      SMIN = max( SMLNUM, EPS*DLANGE( 'M', M, M, A, LDA, DUM ), EPS*DLANGE( 'M', N, N, B, LDB, DUM ) );

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
         for (L = 1; L <= N; L++) { // 60
            if (L < LNEXT) GO TO 60;
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
            for (K = M; K >= 1; K--) { // 50
               if (K > KNEXT) GO TO 50;
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
                  SUML = ddot( M-K1, A( K1, min( K1+1, M ) ), LDA, C( min( K1+1, M ), L1 ), 1 );
                  SUMR = ddot( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 );
                  VEC[1, 1] = C( K1, L1 ) - ( SUML+SGN*SUMR );
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
                  X[1, 1] = ( VEC( 1, 1 )*SCALOC ) / A11;

                  if ( SCALOC != ONE ) {
                     for (J = 1; J <= N; J++) { // 10
                        dscal(M, SCALOC, C( 1, J ), 1 );
                     } // 10
                     SCALE = SCALE*SCALOC;
                  }
                  C[K1, L1] = X( 1, 1 );

               } else if ( L1 == L2 && K1 != K2 ) {

                  SUML = ddot( M-K2, A( K1, min( K2+1, M ) ), LDA, C( min( K2+1, M ), L1 ), 1 );
                  SUMR = ddot( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 );
                  VEC[1, 1] = C( K1, L1 ) - ( SUML+SGN*SUMR );

                  SUML = ddot( M-K2, A( K2, min( K2+1, M ) ), LDA, C( min( K2+1, M ), L1 ), 1 );
                  SUMR = ddot( L1-1, C( K2, 1 ), LDC, B( 1, L1 ), 1 );
                  VEC[2, 1] = C( K2, L1 ) - ( SUML+SGN*SUMR );

                  dlaln2( false , 2, 1, SMIN, ONE, A( K1, K1 ), LDA, ONE, ONE, VEC, 2, -SGN*B( L1, L1 ), ZERO, X, 2, SCALOC, XNORM, IERR );
                  if (IERR != 0) INFO = 1;

                  if ( SCALOC != ONE ) {
                     for (J = 1; J <= N; J++) { // 20
                        dscal(M, SCALOC, C( 1, J ), 1 );
                     } // 20
                     SCALE = SCALE*SCALOC;
                  }
                  C[K1, L1] = X( 1, 1 );
                  C[K2, L1] = X( 2, 1 );

               } else if ( L1 != L2 && K1 == K2 ) {

                  SUML = ddot( M-K1, A( K1, min( K1+1, M ) ), LDA, C( min( K1+1, M ), L1 ), 1 );
                  SUMR = ddot( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 );
                  VEC[1, 1] = SGN*( C( K1, L1 )-( SUML+SGN*SUMR ) );

                  SUML = ddot( M-K1, A( K1, min( K1+1, M ) ), LDA, C( min( K1+1, M ), L2 ), 1 );
                  SUMR = ddot( L1-1, C( K1, 1 ), LDC, B( 1, L2 ), 1 );
                  VEC[2, 1] = SGN*( C( K1, L2 )-( SUML+SGN*SUMR ) );

                  dlaln2( true , 2, 1, SMIN, ONE, B( L1, L1 ), LDB, ONE, ONE, VEC, 2, -SGN*A( K1, K1 ), ZERO, X, 2, SCALOC, XNORM, IERR );
                  if (IERR != 0) INFO = 1;

                  if ( SCALOC != ONE ) {
                     for (J = 1; J <= N; J++) { // 30
                        dscal(M, SCALOC, C( 1, J ), 1 );
                     } // 30
                     SCALE = SCALE*SCALOC;
                  }
                  C[K1, L1] = X( 1, 1 );
                  C[K1, L2] = X( 2, 1 );

               } else if ( L1 != L2 && K1 != K2 ) {

                  SUML = ddot( M-K2, A( K1, min( K2+1, M ) ), LDA, C( min( K2+1, M ), L1 ), 1 );
                  SUMR = ddot( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 );
                  VEC[1, 1] = C( K1, L1 ) - ( SUML+SGN*SUMR );

                  SUML = ddot( M-K2, A( K1, min( K2+1, M ) ), LDA, C( min( K2+1, M ), L2 ), 1 );
                  SUMR = ddot( L1-1, C( K1, 1 ), LDC, B( 1, L2 ), 1 );
                  VEC[1, 2] = C( K1, L2 ) - ( SUML+SGN*SUMR );

                  SUML = ddot( M-K2, A( K2, min( K2+1, M ) ), LDA, C( min( K2+1, M ), L1 ), 1 );
                  SUMR = ddot( L1-1, C( K2, 1 ), LDC, B( 1, L1 ), 1 );
                  VEC[2, 1] = C( K2, L1 ) - ( SUML+SGN*SUMR );

                  SUML = ddot( M-K2, A( K2, min( K2+1, M ) ), LDA, C( min( K2+1, M ), L2 ), 1 );
                  SUMR = ddot( L1-1, C( K2, 1 ), LDC, B( 1, L2 ), 1 );
                  VEC[2, 2] = C( K2, L2 ) - ( SUML+SGN*SUMR );

                  dlasy2( false , false , ISGN, 2, 2, A( K1, K1 ), LDA, B( L1, L1 ), LDB, VEC, 2, SCALOC, X, 2, XNORM, IERR );
                  if (IERR != 0) INFO = 1;

                  if ( SCALOC != ONE ) {
                     for (J = 1; J <= N; J++) { // 40
                        dscal(M, SCALOC, C( 1, J ), 1 );
                     } // 40
                     SCALE = SCALE*SCALOC;
                  }
                  C[K1, L1] = X( 1, 1 );
                  C[K1, L2] = X( 1, 2 );
                  C[K2, L1] = X( 2, 1 );
                  C[K2, L2] = X( 2, 2 );
               }

            } // 50

         } // 60

      } else if ( !NOTRNA && NOTRNB ) {

         // Solve    A**T *X + ISGN*X*B = scale*C.

         // The (K,L)th block of X is determined starting from
         // upper-left corner column by column by

           // A(K,K)**T*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)

         // Where
                    // K-1        T                    L-1
           // R(K,L) = SUM [A(I,K)**T*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)]
                    // I=1                          J=1

         // Start column loop (index = L)
         // L1 (L2): column index of the first (last) row of X(K,L)

         LNEXT = 1;
         for (L = 1; L <= N; L++) { // 120
            if (L < LNEXT) GO TO 120;
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
            for (K = 1; K <= M; K++) { // 110
               if (K < KNEXT) GO TO 110;
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
                  SUML = ddot( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 );
                  SUMR = ddot( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 );
                  VEC[1, 1] = C( K1, L1 ) - ( SUML+SGN*SUMR );
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
                  X[1, 1] = ( VEC( 1, 1 )*SCALOC ) / A11;

                  if ( SCALOC != ONE ) {
                     for (J = 1; J <= N; J++) { // 70
                        dscal(M, SCALOC, C( 1, J ), 1 );
                     } // 70
                     SCALE = SCALE*SCALOC;
                  }
                  C[K1, L1] = X( 1, 1 );

               } else if ( L1 == L2 && K1 != K2 ) {

                  SUML = ddot( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 );
                  SUMR = ddot( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 );
                  VEC[1, 1] = C( K1, L1 ) - ( SUML+SGN*SUMR );

                  SUML = ddot( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 );
                  SUMR = ddot( L1-1, C( K2, 1 ), LDC, B( 1, L1 ), 1 );
                  VEC[2, 1] = C( K2, L1 ) - ( SUML+SGN*SUMR );

                  dlaln2( true , 2, 1, SMIN, ONE, A( K1, K1 ), LDA, ONE, ONE, VEC, 2, -SGN*B( L1, L1 ), ZERO, X, 2, SCALOC, XNORM, IERR );
                  if (IERR != 0) INFO = 1;

                  if ( SCALOC != ONE ) {
                     for (J = 1; J <= N; J++) { // 80
                        dscal(M, SCALOC, C( 1, J ), 1 );
                     } // 80
                     SCALE = SCALE*SCALOC;
                  }
                  C[K1, L1] = X( 1, 1 );
                  C[K2, L1] = X( 2, 1 );

               } else if ( L1 != L2 && K1 == K2 ) {

                  SUML = ddot( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 );
                  SUMR = ddot( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 );
                  VEC[1, 1] = SGN*( C( K1, L1 )-( SUML+SGN*SUMR ) );

                  SUML = ddot( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 );
                  SUMR = ddot( L1-1, C( K1, 1 ), LDC, B( 1, L2 ), 1 );
                  VEC[2, 1] = SGN*( C( K1, L2 )-( SUML+SGN*SUMR ) );

                  dlaln2( true , 2, 1, SMIN, ONE, B( L1, L1 ), LDB, ONE, ONE, VEC, 2, -SGN*A( K1, K1 ), ZERO, X, 2, SCALOC, XNORM, IERR );
                  if (IERR != 0) INFO = 1;

                  if ( SCALOC != ONE ) {
                     for (J = 1; J <= N; J++) { // 90
                        dscal(M, SCALOC, C( 1, J ), 1 );
                     } // 90
                     SCALE = SCALE*SCALOC;
                  }
                  C[K1, L1] = X( 1, 1 );
                  C[K1, L2] = X( 2, 1 );

               } else if ( L1 != L2 && K1 != K2 ) {

                  SUML = ddot( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 );
                  SUMR = ddot( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 );
                  VEC[1, 1] = C( K1, L1 ) - ( SUML+SGN*SUMR );

                  SUML = ddot( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 );
                  SUMR = ddot( L1-1, C( K1, 1 ), LDC, B( 1, L2 ), 1 );
                  VEC[1, 2] = C( K1, L2 ) - ( SUML+SGN*SUMR );

                  SUML = ddot( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 );
                  SUMR = ddot( L1-1, C( K2, 1 ), LDC, B( 1, L1 ), 1 );
                  VEC[2, 1] = C( K2, L1 ) - ( SUML+SGN*SUMR );

                  SUML = ddot( K1-1, A( 1, K2 ), 1, C( 1, L2 ), 1 );
                  SUMR = ddot( L1-1, C( K2, 1 ), LDC, B( 1, L2 ), 1 );
                  VEC[2, 2] = C( K2, L2 ) - ( SUML+SGN*SUMR );

                  dlasy2( true , false , ISGN, 2, 2, A( K1, K1 ), LDA, B( L1, L1 ), LDB, VEC, 2, SCALOC, X, 2, XNORM, IERR );
                  if (IERR != 0) INFO = 1;

                  if ( SCALOC != ONE ) {
                     for (J = 1; J <= N; J++) { // 100
                        dscal(M, SCALOC, C( 1, J ), 1 );
                     } // 100
                     SCALE = SCALE*SCALOC;
                  }
                  C[K1, L1] = X( 1, 1 );
                  C[K1, L2] = X( 1, 2 );
                  C[K2, L1] = X( 2, 1 );
                  C[K2, L2] = X( 2, 2 );
               }

            } // 110
         } // 120

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
         for (L = N; L >= 1; L--) { // 180
            if (L > LNEXT) GO TO 180;
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
            for (K = 1; K <= M; K++) { // 170
               if (K < KNEXT) GO TO 170;
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
                  SUML = ddot( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 );
                  SUMR = ddot( N-L1, C( K1, min( L1+1, N ) ), LDC, B( L1, min( L1+1, N ) ), LDB );
                  VEC[1, 1] = C( K1, L1 ) - ( SUML+SGN*SUMR );
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
                  X[1, 1] = ( VEC( 1, 1 )*SCALOC ) / A11;

                  if ( SCALOC != ONE ) {
                     for (J = 1; J <= N; J++) { // 130
                        dscal(M, SCALOC, C( 1, J ), 1 );
                     } // 130
                     SCALE = SCALE*SCALOC;
                  }
                  C[K1, L1] = X( 1, 1 );

               } else if ( L1 == L2 && K1 != K2 ) {

                  SUML = ddot( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 );
                  SUMR = ddot( N-L2, C( K1, min( L2+1, N ) ), LDC, B( L1, min( L2+1, N ) ), LDB );
                  VEC[1, 1] = C( K1, L1 ) - ( SUML+SGN*SUMR );

                  SUML = ddot( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 );
                  SUMR = ddot( N-L2, C( K2, min( L2+1, N ) ), LDC, B( L1, min( L2+1, N ) ), LDB );
                  VEC[2, 1] = C( K2, L1 ) - ( SUML+SGN*SUMR );

                  dlaln2( true , 2, 1, SMIN, ONE, A( K1, K1 ), LDA, ONE, ONE, VEC, 2, -SGN*B( L1, L1 ), ZERO, X, 2, SCALOC, XNORM, IERR );
                  if (IERR != 0) INFO = 1;

                  if ( SCALOC != ONE ) {
                     for (J = 1; J <= N; J++) { // 140
                        dscal(M, SCALOC, C( 1, J ), 1 );
                     } // 140
                     SCALE = SCALE*SCALOC;
                  }
                  C[K1, L1] = X( 1, 1 );
                  C[K2, L1] = X( 2, 1 );

               } else if ( L1 != L2 && K1 == K2 ) {

                  SUML = ddot( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 );
                  SUMR = ddot( N-L2, C( K1, min( L2+1, N ) ), LDC, B( L1, min( L2+1, N ) ), LDB );
                  VEC[1, 1] = SGN*( C( K1, L1 )-( SUML+SGN*SUMR ) );

                  SUML = ddot( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 );
                  SUMR = ddot( N-L2, C( K1, min( L2+1, N ) ), LDC, B( L2, min( L2+1, N ) ), LDB );
                  VEC[2, 1] = SGN*( C( K1, L2 )-( SUML+SGN*SUMR ) );

                  dlaln2( false , 2, 1, SMIN, ONE, B( L1, L1 ), LDB, ONE, ONE, VEC, 2, -SGN*A( K1, K1 ), ZERO, X, 2, SCALOC, XNORM, IERR );
                  if (IERR != 0) INFO = 1;

                  if ( SCALOC != ONE ) {
                     for (J = 1; J <= N; J++) { // 150
                        dscal(M, SCALOC, C( 1, J ), 1 );
                     } // 150
                     SCALE = SCALE*SCALOC;
                  }
                  C[K1, L1] = X( 1, 1 );
                  C[K1, L2] = X( 2, 1 );

               } else if ( L1 != L2 && K1 != K2 ) {

                  SUML = ddot( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 );
                  SUMR = ddot( N-L2, C( K1, min( L2+1, N ) ), LDC, B( L1, min( L2+1, N ) ), LDB );
                  VEC[1, 1] = C( K1, L1 ) - ( SUML+SGN*SUMR );

                  SUML = ddot( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 );
                  SUMR = ddot( N-L2, C( K1, min( L2+1, N ) ), LDC, B( L2, min( L2+1, N ) ), LDB );
                  VEC[1, 2] = C( K1, L2 ) - ( SUML+SGN*SUMR );

                  SUML = ddot( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 );
                  SUMR = ddot( N-L2, C( K2, min( L2+1, N ) ), LDC, B( L1, min( L2+1, N ) ), LDB );
                  VEC[2, 1] = C( K2, L1 ) - ( SUML+SGN*SUMR );

                  SUML = ddot( K1-1, A( 1, K2 ), 1, C( 1, L2 ), 1 );
                  SUMR = ddot( N-L2, C( K2, min( L2+1, N ) ), LDC, B( L2, min( L2+1, N ) ), LDB );
                  VEC[2, 2] = C( K2, L2 ) - ( SUML+SGN*SUMR );

                  dlasy2( true , true , ISGN, 2, 2, A( K1, K1 ), LDA, B( L1, L1 ), LDB, VEC, 2, SCALOC, X, 2, XNORM, IERR );
                  if (IERR != 0) INFO = 1;

                  if ( SCALOC != ONE ) {
                     for (J = 1; J <= N; J++) { // 160
                        dscal(M, SCALOC, C( 1, J ), 1 );
                     } // 160
                     SCALE = SCALE*SCALOC;
                  }
                  C[K1, L1] = X( 1, 1 );
                  C[K1, L2] = X( 1, 2 );
                  C[K2, L1] = X( 2, 1 );
                  C[K2, L2] = X( 2, 2 );
               }

            } // 170
         } // 180

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
         for (L = N; L >= 1; L--) { // 240
            if (L > LNEXT) GO TO 240;
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
            for (K = M; K >= 1; K--) { // 230
               if (K > KNEXT) GO TO 230;
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
                  SUML = ddot( M-K1, A( K1, min( K1+1, M ) ), LDA, C( min( K1+1, M ), L1 ), 1 )                   SUMR = ddot( N-L1, C( K1, min( L1+1, N ) ), LDC, B( L1, min( L1+1, N ) ), LDB );
                  VEC[1, 1] = C( K1, L1 ) - ( SUML+SGN*SUMR );
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
                  X[1, 1] = ( VEC( 1, 1 )*SCALOC ) / A11;

                  if ( SCALOC != ONE ) {
                     for (J = 1; J <= N; J++) { // 190
                        dscal(M, SCALOC, C( 1, J ), 1 );
                     } // 190
                     SCALE = SCALE*SCALOC;
                  }
                  C[K1, L1] = X( 1, 1 );

               } else if ( L1 == L2 && K1 != K2 ) {

                  SUML = ddot( M-K2, A( K1, min( K2+1, M ) ), LDA, C( min( K2+1, M ), L1 ), 1 )                   SUMR = ddot( N-L2, C( K1, min( L2+1, N ) ), LDC, B( L1, min( L2+1, N ) ), LDB );
                  VEC[1, 1] = C( K1, L1 ) - ( SUML+SGN*SUMR );

                  SUML = ddot( M-K2, A( K2, min( K2+1, M ) ), LDA, C( min( K2+1, M ), L1 ), 1 )                   SUMR = ddot( N-L2, C( K2, min( L2+1, N ) ), LDC, B( L1, min( L2+1, N ) ), LDB );
                  VEC[2, 1] = C( K2, L1 ) - ( SUML+SGN*SUMR );

                  dlaln2( false , 2, 1, SMIN, ONE, A( K1, K1 ), LDA, ONE, ONE, VEC, 2, -SGN*B( L1, L1 ), ZERO, X, 2, SCALOC, XNORM, IERR );
                  if (IERR != 0) INFO = 1;

                  if ( SCALOC != ONE ) {
                     for (J = 1; J <= N; J++) { // 200
                        dscal(M, SCALOC, C( 1, J ), 1 );
                     } // 200
                     SCALE = SCALE*SCALOC;
                  }
                  C[K1, L1] = X( 1, 1 );
                  C[K2, L1] = X( 2, 1 );

               } else if ( L1 != L2 && K1 == K2 ) {

                  SUML = ddot( M-K1, A( K1, min( K1+1, M ) ), LDA, C( min( K1+1, M ), L1 ), 1 )                   SUMR = ddot( N-L2, C( K1, min( L2+1, N ) ), LDC, B( L1, min( L2+1, N ) ), LDB );
                  VEC[1, 1] = SGN*( C( K1, L1 )-( SUML+SGN*SUMR ) );

                  SUML = ddot( M-K1, A( K1, min( K1+1, M ) ), LDA, C( min( K1+1, M ), L2 ), 1 )                   SUMR = ddot( N-L2, C( K1, min( L2+1, N ) ), LDC, B( L2, min( L2+1, N ) ), LDB );
                  VEC[2, 1] = SGN*( C( K1, L2 )-( SUML+SGN*SUMR ) );

                  dlaln2( false , 2, 1, SMIN, ONE, B( L1, L1 ), LDB, ONE, ONE, VEC, 2, -SGN*A( K1, K1 ), ZERO, X, 2, SCALOC, XNORM, IERR );
                  if (IERR != 0) INFO = 1;

                  if ( SCALOC != ONE ) {
                     for (J = 1; J <= N; J++) { // 210
                        dscal(M, SCALOC, C( 1, J ), 1 );
                     } // 210
                     SCALE = SCALE*SCALOC;
                  }
                  C[K1, L1] = X( 1, 1 );
                  C[K1, L2] = X( 2, 1 );

               } else if ( L1 != L2 && K1 != K2 ) {

                  SUML = ddot( M-K2, A( K1, min( K2+1, M ) ), LDA, C( min( K2+1, M ), L1 ), 1 )                   SUMR = ddot( N-L2, C( K1, min( L2+1, N ) ), LDC, B( L1, min( L2+1, N ) ), LDB );
                  VEC[1, 1] = C( K1, L1 ) - ( SUML+SGN*SUMR );

                  SUML = ddot( M-K2, A( K1, min( K2+1, M ) ), LDA, C( min( K2+1, M ), L2 ), 1 )                   SUMR = ddot( N-L2, C( K1, min( L2+1, N ) ), LDC, B( L2, min( L2+1, N ) ), LDB );
                  VEC[1, 2] = C( K1, L2 ) - ( SUML+SGN*SUMR );

                  SUML = ddot( M-K2, A( K2, min( K2+1, M ) ), LDA, C( min( K2+1, M ), L1 ), 1 )                   SUMR = ddot( N-L2, C( K2, min( L2+1, N ) ), LDC, B( L1, min( L2+1, N ) ), LDB );
                  VEC[2, 1] = C( K2, L1 ) - ( SUML+SGN*SUMR );

                  SUML = ddot( M-K2, A( K2, min( K2+1, M ) ), LDA, C( min( K2+1, M ), L2 ), 1 )                   SUMR = ddot( N-L2, C( K2, min( L2+1, N ) ), LDC, B( L2, min( L2+1, N ) ), LDB );
                  VEC[2, 2] = C( K2, L2 ) - ( SUML+SGN*SUMR );

                  dlasy2( false , true , ISGN, 2, 2, A( K1, K1 ), LDA, B( L1, L1 ), LDB, VEC, 2, SCALOC, X, 2, XNORM, IERR );
                  if (IERR != 0) INFO = 1;

                  if ( SCALOC != ONE ) {
                     for (J = 1; J <= N; J++) { // 220
                        dscal(M, SCALOC, C( 1, J ), 1 );
                     } // 220
                     SCALE = SCALE*SCALOC;
                  }
                  C[K1, L1] = X( 1, 1 );
                  C[K1, L2] = X( 1, 2 );
                  C[K2, L1] = X( 2, 1 );
                  C[K2, L2] = X( 2, 2 );
               }

            } // 230
         } // 240

      }

      return;
      }
