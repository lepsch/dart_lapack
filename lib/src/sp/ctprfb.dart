      void ctprfb(final int SIDE, final int TRANS, final int DIRECT, final int STOREV, final int M, final int N, final int K, final int L, final Matrix<double> V, final int LDV, final Matrix<double> T, final int LDT, final Matrix<double> A, final int LDA, final Matrix<double> B, final int LDB, final Array<double> _WORK, final int LDWORK,) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String    DIRECT, SIDE, STOREV, TRANS;
      int       K, L, LDA, LDB, LDT, LDV, LDWORK, M, N;
      Complex   A( LDA, * ), B( LDB, * ), T( LDT, * ), V( LDV, * ), WORK( LDWORK, * );
      // ..

// ==========================================================================

      // .. Parameters ..
      Complex   ONE, ZERO;
      const     ONE = (1.0,0.0), ZERO = (0.0,0.0) ;
      int       I, J, MP, NP, KP;
      bool      LEFT, FORWARD, COLUMN, RIGHT, BACKWARD, ROW;
      // ..
      // .. External Functions ..
      //- bool      lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM, CTRMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG

      // Quick return if possible

      if (M <= 0 || N <= 0 || K <= 0 || L < 0) return;

      if ( lsame( STOREV, 'C' ) ) {
         COLUMN = true;
         ROW = false;
      } else if ( lsame( STOREV, 'R' ) ) {
         COLUMN = false;
         ROW = true;
      } else {
         COLUMN = false;
         ROW = false;
      }

      if ( lsame( SIDE, 'L' ) ) {
         LEFT = true;
         RIGHT = false;
      } else if ( lsame( SIDE, 'R' ) ) {
         LEFT = false;
         RIGHT = true;
      } else {
         LEFT = false;
         RIGHT = false;
      }

      if ( lsame( DIRECT, 'F' ) ) {
         FORWARD = true;
         BACKWARD = false;
      } else if ( lsame( DIRECT, 'B' ) ) {
         FORWARD = false;
         BACKWARD = true;
      } else {
         FORWARD = false;
         BACKWARD = false;
      }

// ---------------------------------------------------------------------------

      if ( COLUMN && FORWARD && LEFT  ) {

// ---------------------------------------------------------------------------

         // Let  W =  [ I ]    (K-by-K)
         //           [ V ]    (M-by-K)

         // Form  H C  or  H**H C  where  C = [ A ]  (K-by-N)
         //                                   [ B ]  (M-by-N)

         // H = I - W T W**H          or  H**H = I - W T**H W**H

         // A = A -   T (A + V**H B)  or  A = A -   T**H (A + V**H B)
         // B = B - V T (A + V**H B)  or  B = B - V T**H (A + V**H B)

// ---------------------------------------------------------------------------

         MP = min( M-L+1, M );
         KP = min( L+1, K );

         for (J = 1; J <= N; J++) {
            for (I = 1; I <= L; I++) {
               WORK[I][J] = B( M-L+I, J );
            }
         }
         ctrmm('L', 'U', 'C', 'N', L, N, ONE, V( MP, 1 ), LDV, WORK, LDWORK );
         cgemm('C', 'N', L, N, M-L, ONE, V, LDV, B, LDB, ONE, WORK, LDWORK );
         cgemm('C', 'N', K-L, N, M, ONE, V( 1, KP ), LDV, B, LDB, ZERO, WORK( KP, 1 ), LDWORK );

         for (J = 1; J <= N; J++) {
            for (I = 1; I <= K; I++) {
               WORK[I][J] = WORK( I, J ) + A( I, J );
            }
         }

         ctrmm('L', 'U', TRANS, 'N', K, N, ONE, T, LDT, WORK, LDWORK );

         for (J = 1; J <= N; J++) {
            for (I = 1; I <= K; I++) {
               A[I][J] = A( I, J ) - WORK( I, J );
            }
         }

         cgemm('N', 'N', M-L, N, K, -ONE, V, LDV, WORK, LDWORK, ONE, B, LDB );
         cgemm('N', 'N', L, N, K-L, -ONE, V( MP, KP ), LDV, WORK( KP, 1 ), LDWORK, ONE, B( MP, 1 ),  LDB );
         ctrmm('L', 'U', 'N', 'N', L, N, ONE, V( MP, 1 ), LDV, WORK, LDWORK );
         for (J = 1; J <= N; J++) {
            for (I = 1; I <= L; I++) {
               B[M-L+I][J] = B( M-L+I, J ) - WORK( I, J );
            }
         }

// ---------------------------------------------------------------------------

      } else if ( COLUMN && FORWARD && RIGHT ) {

// ---------------------------------------------------------------------------

         // Let  W =  [ I ]    (K-by-K)
         //           [ V ]    (N-by-K)

         // Form  C H or  C H**H  where  C = [ A B ] (A is M-by-K, B is M-by-N)

         // H = I - W T W**H          or  H**H = I - W T**H W**H

         // A = A - (A + B V) T      or  A = A - (A + B V) T**H
         // B = B - (A + B V) T V**H  or  B = B - (A + B V) T**H V**H

// ---------------------------------------------------------------------------

         NP = min( N-L+1, N );
         KP = min( L+1, K );

         for (J = 1; J <= L; J++) {
            for (I = 1; I <= M; I++) {
               WORK[I][J] = B( I, N-L+J );
            }
         }
         ctrmm('R', 'U', 'N', 'N', M, L, ONE, V( NP, 1 ), LDV, WORK, LDWORK );
         cgemm('N', 'N', M, L, N-L, ONE, B, LDB, V, LDV, ONE, WORK, LDWORK );
         cgemm('N', 'N', M, K-L, N, ONE, B, LDB, V( 1, KP ), LDV, ZERO, WORK( 1, KP ), LDWORK );

         for (J = 1; J <= K; J++) {
            for (I = 1; I <= M; I++) {
               WORK[I][J] = WORK( I, J ) + A( I, J );
            }
         }

         ctrmm('R', 'U', TRANS, 'N', M, K, ONE, T, LDT, WORK, LDWORK );

         for (J = 1; J <= K; J++) {
            for (I = 1; I <= M; I++) {
               A[I][J] = A( I, J ) - WORK( I, J );
            }
         }

         cgemm('N', 'C', M, N-L, K, -ONE, WORK, LDWORK, V, LDV, ONE, B, LDB );
         cgemm('N', 'C', M, L, K-L, -ONE, WORK( 1, KP ), LDWORK, V( NP, KP ), LDV, ONE, B( 1, NP ), LDB );
         ctrmm('R', 'U', 'C', 'N', M, L, ONE, V( NP, 1 ), LDV, WORK, LDWORK );
         for (J = 1; J <= L; J++) {
            for (I = 1; I <= M; I++) {
               B[I][N-L+J] = B( I, N-L+J ) - WORK( I, J );
            }
         }

// ---------------------------------------------------------------------------

      } else if ( COLUMN && BACKWARD && LEFT ) {

// ---------------------------------------------------------------------------

         // Let  W =  [ V ]    (M-by-K)
         //           [ I ]    (K-by-K)

         // Form  H C  or  H**H C  where  C = [ B ]  (M-by-N)
         //                                   [ A ]  (K-by-N)

         // H = I - W T W**H          or  H**H = I - W T**H W**H

         // A = A -   T (A + V**H B)  or  A = A -   T**H (A + V**H B)
         // B = B - V T (A + V**H B)  or  B = B - V T**H (A + V**H B)

// ---------------------------------------------------------------------------

         MP = min( L+1, M );
         KP = min( K-L+1, K );

         for (J = 1; J <= N; J++) {
            for (I = 1; I <= L; I++) {
               WORK[K-L+I][J] = B( I, J );
            }
         }

         ctrmm('L', 'L', 'C', 'N', L, N, ONE, V( 1, KP ), LDV, WORK( KP, 1 ), LDWORK );
         cgemm('C', 'N', L, N, M-L, ONE, V( MP, KP ), LDV, B( MP, 1 ), LDB, ONE, WORK( KP, 1 ), LDWORK );
         cgemm('C', 'N', K-L, N, M, ONE, V, LDV, B, LDB, ZERO, WORK, LDWORK );

         for (J = 1; J <= N; J++) {
            for (I = 1; I <= K; I++) {
               WORK[I][J] = WORK( I, J ) + A( I, J );
            }
         }

         ctrmm('L', 'L', TRANS, 'N', K, N, ONE, T, LDT, WORK, LDWORK );

         for (J = 1; J <= N; J++) {
            for (I = 1; I <= K; I++) {
               A[I][J] = A( I, J ) - WORK( I, J );
            }
         }

         cgemm('N', 'N', M-L, N, K, -ONE, V( MP, 1 ), LDV, WORK, LDWORK, ONE, B( MP, 1 ), LDB );
         cgemm('N', 'N', L, N, K-L, -ONE, V, LDV, WORK, LDWORK, ONE, B,  LDB );
         ctrmm('L', 'L', 'N', 'N', L, N, ONE, V( 1, KP ), LDV, WORK( KP, 1 ), LDWORK );
         for (J = 1; J <= N; J++) {
            for (I = 1; I <= L; I++) {
               B[I][J] = B( I, J ) - WORK( K-L+I, J );
            }
         }

// ---------------------------------------------------------------------------

      } else if ( COLUMN && BACKWARD && RIGHT ) {

// ---------------------------------------------------------------------------

         // Let  W =  [ V ]    (N-by-K)
         //           [ I ]    (K-by-K)

         // Form  C H  or  C H**H  where  C = [ B A ] (B is M-by-N, A is M-by-K)

         // H = I - W T W**H          or  H**H = I - W T**H W**H

         // A = A - (A + B V) T      or  A = A - (A + B V) T**H
         // B = B - (A + B V) T V**H  or  B = B - (A + B V) T**H V**H

// ---------------------------------------------------------------------------

         NP = min( L+1, N );
         KP = min( K-L+1, K );

         for (J = 1; J <= L; J++) {
            for (I = 1; I <= M; I++) {
               WORK[I][K-L+J] = B( I, J );
            }
         }
         ctrmm('R', 'L', 'N', 'N', M, L, ONE, V( 1, KP ), LDV, WORK( 1, KP ), LDWORK );
         cgemm('N', 'N', M, L, N-L, ONE, B( 1, NP ), LDB, V( NP, KP ), LDV, ONE, WORK( 1, KP ), LDWORK );
         cgemm('N', 'N', M, K-L, N, ONE, B, LDB, V, LDV, ZERO, WORK, LDWORK );

         for (J = 1; J <= K; J++) {
            for (I = 1; I <= M; I++) {
               WORK[I][J] = WORK( I, J ) + A( I, J );
            }
         }

         ctrmm('R', 'L', TRANS, 'N', M, K, ONE, T, LDT, WORK, LDWORK );

         for (J = 1; J <= K; J++) {
            for (I = 1; I <= M; I++) {
               A[I][J] = A( I, J ) - WORK( I, J );
            }
         }

         cgemm('N', 'C', M, N-L, K, -ONE, WORK, LDWORK, V( NP, 1 ), LDV, ONE, B( 1, NP ), LDB );
         cgemm('N', 'C', M, L, K-L, -ONE, WORK, LDWORK, V, LDV, ONE, B, LDB );
         ctrmm('R', 'L', 'C', 'N', M, L, ONE, V( 1, KP ), LDV, WORK( 1, KP ), LDWORK );
         for (J = 1; J <= L; J++) {
            for (I = 1; I <= M; I++) {
               B[I][J] = B( I, J ) - WORK( I, K-L+J );
            }
         }

// ---------------------------------------------------------------------------

      } else if ( ROW && FORWARD && LEFT ) {

// ---------------------------------------------------------------------------

         // Let  W =  [ I V ] ( I is K-by-K, V is K-by-M )

         // Form  H C  or  H**H C  where  C = [ A ]  (K-by-N)
         //                                   [ B ]  (M-by-N)

         // H = I - W**H T W          or  H**H = I - W**H T**H W

         // A = A -     T (A + V B)  or  A = A -     T**H (A + V B)
         // B = B - V**H T (A + V B)  or  B = B - V**H T**H (A + V B)

// ---------------------------------------------------------------------------

         MP = min( M-L+1, M );
         KP = min( L+1, K );

         for (J = 1; J <= N; J++) {
            for (I = 1; I <= L; I++) {
               WORK[I][J] = B( M-L+I, J );
            }
         }
         ctrmm('L', 'L', 'N', 'N', L, N, ONE, V( 1, MP ), LDV, WORK, LDB );
         cgemm('N', 'N', L, N, M-L, ONE, V, LDV,B, LDB, ONE, WORK, LDWORK );
         cgemm('N', 'N', K-L, N, M, ONE, V( KP, 1 ), LDV, B, LDB, ZERO, WORK( KP, 1 ), LDWORK );

         for (J = 1; J <= N; J++) {
            for (I = 1; I <= K; I++) {
               WORK[I][J] = WORK( I, J ) + A( I, J );
            }
         }

         ctrmm('L', 'U', TRANS, 'N', K, N, ONE, T, LDT, WORK, LDWORK );

         for (J = 1; J <= N; J++) {
            for (I = 1; I <= K; I++) {
               A[I][J] = A( I, J ) - WORK( I, J );
            }
         }

         cgemm('C', 'N', M-L, N, K, -ONE, V, LDV, WORK, LDWORK, ONE, B, LDB );
         cgemm('C', 'N', L, N, K-L, -ONE, V( KP, MP ), LDV, WORK( KP, 1 ), LDWORK, ONE, B( MP, 1 ), LDB );
         ctrmm('L', 'L', 'C', 'N', L, N, ONE, V( 1, MP ), LDV, WORK, LDWORK );
         for (J = 1; J <= N; J++) {
            for (I = 1; I <= L; I++) {
               B[M-L+I][J] = B( M-L+I, J ) - WORK( I, J );
            }
         }

// ---------------------------------------------------------------------------

      } else if ( ROW && FORWARD && RIGHT ) {

// ---------------------------------------------------------------------------

         // Let  W =  [ I V ] ( I is K-by-K, V is K-by-N )

         // Form  C H  or  C H**H  where  C = [ A B ] (A is M-by-K, B is M-by-N)

         // H = I - W**H T W            or  H**H = I - W**H T**H W

         // A = A - (A + B V**H) T      or  A = A - (A + B V**H) T**H
         // B = B - (A + B V**H) T V    or  B = B - (A + B V**H) T**H V

// ---------------------------------------------------------------------------

         NP = min( N-L+1, N );
         KP = min( L+1, K );

         for (J = 1; J <= L; J++) {
            for (I = 1; I <= M; I++) {
               WORK[I][J] = B( I, N-L+J );
            }
         }
         ctrmm('R', 'L', 'C', 'N', M, L, ONE, V( 1, NP ), LDV, WORK, LDWORK );
         cgemm('N', 'C', M, L, N-L, ONE, B, LDB, V, LDV, ONE, WORK, LDWORK );
         cgemm('N', 'C', M, K-L, N, ONE, B, LDB, V( KP, 1 ), LDV, ZERO, WORK( 1, KP ), LDWORK );

         for (J = 1; J <= K; J++) {
            for (I = 1; I <= M; I++) {
               WORK[I][J] = WORK( I, J ) + A( I, J );
            }
         }

         ctrmm('R', 'U', TRANS, 'N', M, K, ONE, T, LDT, WORK, LDWORK );

         for (J = 1; J <= K; J++) {
            for (I = 1; I <= M; I++) {
               A[I][J] = A( I, J ) - WORK( I, J );
            }
         }

         cgemm('N', 'N', M, N-L, K, -ONE, WORK, LDWORK, V, LDV, ONE, B, LDB );
         cgemm('N', 'N', M, L, K-L, -ONE, WORK( 1, KP ), LDWORK, V( KP, NP ), LDV, ONE, B( 1, NP ), LDB );
         ctrmm('R', 'L', 'N', 'N', M, L, ONE, V( 1, NP ), LDV, WORK, LDWORK );
         for (J = 1; J <= L; J++) {
            for (I = 1; I <= M; I++) {
               B[I][N-L+J] = B( I, N-L+J ) - WORK( I, J );
            }
         }

// ---------------------------------------------------------------------------

      } else if ( ROW && BACKWARD && LEFT ) {

// ---------------------------------------------------------------------------

         // Let  W =  [ V I ] ( I is K-by-K, V is K-by-M )

         // Form  H C  or  H**H C  where  C = [ B ]  (M-by-N)
         //                                   [ A ]  (K-by-N)

         // H = I - W**H T W          or  H**H = I - W**H T**H W

         // A = A -     T (A + V B)  or  A = A -     T**H (A + V B)
         // B = B - V**H T (A + V B)  or  B = B - V**H T**H (A + V B)

// ---------------------------------------------------------------------------

         MP = min( L+1, M );
         KP = min( K-L+1, K );

         for (J = 1; J <= N; J++) {
            for (I = 1; I <= L; I++) {
               WORK[K-L+I][J] = B( I, J );
            }
         }
         ctrmm('L', 'U', 'N', 'N', L, N, ONE, V( KP, 1 ), LDV, WORK( KP, 1 ), LDWORK );
         cgemm('N', 'N', L, N, M-L, ONE, V( KP, MP ), LDV, B( MP, 1 ), LDB, ONE, WORK( KP, 1 ), LDWORK );
         cgemm('N', 'N', K-L, N, M, ONE, V, LDV, B, LDB, ZERO, WORK, LDWORK );

         for (J = 1; J <= N; J++) {
            for (I = 1; I <= K; I++) {
               WORK[I][J] = WORK( I, J ) + A( I, J );
            }
         }

         ctrmm('L', 'L ', TRANS, 'N', K, N, ONE, T, LDT, WORK, LDWORK );

         for (J = 1; J <= N; J++) {
            for (I = 1; I <= K; I++) {
               A[I][J] = A( I, J ) - WORK( I, J );
            }
         }

         cgemm('C', 'N', M-L, N, K, -ONE, V( 1, MP ), LDV, WORK, LDWORK, ONE, B( MP, 1 ), LDB );
         cgemm('C', 'N', L, N, K-L, -ONE, V, LDV, WORK, LDWORK, ONE, B, LDB );
         ctrmm('L', 'U', 'C', 'N', L, N, ONE, V( KP, 1 ), LDV, WORK( KP, 1 ), LDWORK );
         for (J = 1; J <= N; J++) {
            for (I = 1; I <= L; I++) {
               B[I][J] = B( I, J ) - WORK( K-L+I, J );
            }
         }

// ---------------------------------------------------------------------------

      } else if ( ROW && BACKWARD && RIGHT ) {

// ---------------------------------------------------------------------------

         // Let  W =  [ V I ] ( I is K-by-K, V is K-by-N )

         // Form  C H  or  C H**H  where  C = [ B A ] (A is M-by-K, B is M-by-N)

         // H = I - W**H T W            or  H**H = I - W**H T**H W

         // A = A - (A + B V**H) T      or  A = A - (A + B V**H) T**H
         // B = B - (A + B V**H) T V    or  B = B - (A + B V**H) T**H V

// ---------------------------------------------------------------------------

         NP = min( L+1, N );
         KP = min( K-L+1, K );

         for (J = 1; J <= L; J++) {
            for (I = 1; I <= M; I++) {
               WORK[I][K-L+J] = B( I, J );
            }
         }
         ctrmm('R', 'U', 'C', 'N', M, L, ONE, V( KP, 1 ), LDV, WORK( 1, KP ), LDWORK );
         cgemm('N', 'C', M, L, N-L, ONE, B( 1, NP ), LDB, V( KP, NP ), LDV, ONE, WORK( 1, KP ), LDWORK );
         cgemm('N', 'C', M, K-L, N, ONE, B, LDB, V, LDV, ZERO, WORK, LDWORK );

         for (J = 1; J <= K; J++) {
            for (I = 1; I <= M; I++) {
               WORK[I][J] = WORK( I, J ) + A( I, J );
            }
         }

         ctrmm('R', 'L', TRANS, 'N', M, K, ONE, T, LDT, WORK, LDWORK );

         for (J = 1; J <= K; J++) {
            for (I = 1; I <= M; I++) {
               A[I][J] = A( I, J ) - WORK( I, J );
            }
         }

         cgemm('N', 'N', M, N-L, K, -ONE, WORK, LDWORK, V( 1, NP ), LDV, ONE, B( 1, NP ), LDB );
         cgemm('N', 'N', M, L, K-L , -ONE, WORK, LDWORK, V, LDV, ONE, B, LDB );
         ctrmm('R', 'U', 'N', 'N', M, L, ONE, V( KP, 1 ), LDV, WORK( 1, KP ), LDWORK );
         for (J = 1; J <= L; J++) {
            for (I = 1; I <= M; I++) {
               B[I][J] = B( I, J ) - WORK( I, K-L+J );
            }
         }

      }

      }
