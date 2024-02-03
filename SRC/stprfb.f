      void stprfb(SIDE, TRANS, DIRECT, STOREV, M, N, K, L, V, LDV, T, LDT, A, LDA, B, LDB, WORK, LDWORK ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String    DIRECT, SIDE, STOREV, TRANS;
      int       K, L, LDA, LDB, LDT, LDV, LDWORK, M, N;
      // ..
      // .. Array Arguments ..
      REAL   A( LDA, * ), B( LDB, * ), T( LDT, * ), V( LDV, * ), WORK( LDWORK, * );
      // ..

// ==========================================================================

      // .. Parameters ..
      REAL   ONE, ZERO;
      const     ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int       I, J, MP, NP, KP;
      bool      LEFT, FORWARD, COLUMN, RIGHT, BACKWARD, ROW;
      // ..
      // .. External Functions ..
      //- bool      LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, STRMM
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if (M <= 0 || N <= 0 || K <= 0 || L < 0) return;

      if ( LSAME( STOREV, 'C' ) ) {
         COLUMN = true;
         ROW = false;
      } else if ( LSAME( STOREV, 'R' ) ) {
         COLUMN = false;
         ROW = true;
      } else {
         COLUMN = false;
         ROW = false;
      }

      if ( LSAME( SIDE, 'L' ) ) {
         LEFT = true;
         RIGHT = false;
      } else if ( LSAME( SIDE, 'R' ) ) {
         LEFT = false;
         RIGHT = true;
      } else {
         LEFT = false;
         RIGHT = false;
      }

      if ( LSAME( DIRECT, 'F' ) ) {
         FORWARD = true;
         BACKWARD = false;
      } else if ( LSAME( DIRECT, 'B' ) ) {
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
                   // [ V ]    (M-by-K)

         // Form  H C  or  H**T C  where  C = [ A ]  (K-by-N)
                                           // [ B ]  (M-by-N)

         // H = I - W T W**T          or  H**T = I - W T**T W**T

         // A = A -   T (A + V**T B)  or  A = A -   T**T (A + V**T B)
         // B = B - V T (A + V**T B)  or  B = B - V T**T (A + V**T B)

// ---------------------------------------------------------------------------

         MP = min( M-L+1, M );
         KP = min( L+1, K );

         for (J = 1; J <= N; J++) {
            for (I = 1; I <= L; I++) {
               WORK( I, J ) = B( M-L+I, J );
            }
         }
         strmm('L', 'U', 'T', 'N', L, N, ONE, V( MP, 1 ), LDV, WORK, LDWORK );
         sgemm('T', 'N', L, N, M-L, ONE, V, LDV, B, LDB, ONE, WORK, LDWORK );
         sgemm('T', 'N', K-L, N, M, ONE, V( 1, KP ), LDV, B, LDB, ZERO, WORK( KP, 1 ), LDWORK );

         for (J = 1; J <= N; J++) {
            for (I = 1; I <= K; I++) {
               WORK( I, J ) = WORK( I, J ) + A( I, J );
            }
         }

         strmm('L', 'U', TRANS, 'N', K, N, ONE, T, LDT, WORK, LDWORK );

         for (J = 1; J <= N; J++) {
            for (I = 1; I <= K; I++) {
               A( I, J ) = A( I, J ) - WORK( I, J );
            }
         }

         sgemm('N', 'N', M-L, N, K, -ONE, V, LDV, WORK, LDWORK, ONE, B, LDB );
         sgemm('N', 'N', L, N, K-L, -ONE, V( MP, KP ), LDV, WORK( KP, 1 ), LDWORK, ONE, B( MP, 1 ),  LDB );
         strmm('L', 'U', 'N', 'N', L, N, ONE, V( MP, 1 ), LDV, WORK, LDWORK );
         for (J = 1; J <= N; J++) {
            for (I = 1; I <= L; I++) {
               B( M-L+I, J ) = B( M-L+I, J ) - WORK( I, J );
            }
         }

// ---------------------------------------------------------------------------

      } else if ( COLUMN && FORWARD && RIGHT ) {

// ---------------------------------------------------------------------------

         // Let  W =  [ I ]    (K-by-K)
                   // [ V ]    (N-by-K)

         // Form  C H or  C H**T  where  C = [ A B ] (A is M-by-K, B is M-by-N)

         // H = I - W T W**T          or  H**T = I - W T**T W**T

         // A = A - (A + B V) T       or  A = A - (A + B V) T**T
         // B = B - (A + B V) T V**T  or  B = B - (A + B V) T**T V**T

// ---------------------------------------------------------------------------

         NP = min( N-L+1, N );
         KP = min( L+1, K );

         for (J = 1; J <= L; J++) {
            for (I = 1; I <= M; I++) {
               WORK( I, J ) = B( I, N-L+J );
            }
         }
         strmm('R', 'U', 'N', 'N', M, L, ONE, V( NP, 1 ), LDV, WORK, LDWORK );
         sgemm('N', 'N', M, L, N-L, ONE, B, LDB, V, LDV, ONE, WORK, LDWORK );
         sgemm('N', 'N', M, K-L, N, ONE, B, LDB, V( 1, KP ), LDV, ZERO, WORK( 1, KP ), LDWORK );

         for (J = 1; J <= K; J++) {
            for (I = 1; I <= M; I++) {
               WORK( I, J ) = WORK( I, J ) + A( I, J );
            }
         }

         strmm('R', 'U', TRANS, 'N', M, K, ONE, T, LDT, WORK, LDWORK );

         for (J = 1; J <= K; J++) {
            for (I = 1; I <= M; I++) {
               A( I, J ) = A( I, J ) - WORK( I, J );
            }
         }

         sgemm('N', 'T', M, N-L, K, -ONE, WORK, LDWORK, V, LDV, ONE, B, LDB );
         sgemm('N', 'T', M, L, K-L, -ONE, WORK( 1, KP ), LDWORK, V( NP, KP ), LDV, ONE, B( 1, NP ), LDB );
         strmm('R', 'U', 'T', 'N', M, L, ONE, V( NP, 1 ), LDV, WORK, LDWORK );
         for (J = 1; J <= L; J++) {
            for (I = 1; I <= M; I++) {
               B( I, N-L+J ) = B( I, N-L+J ) - WORK( I, J );
            }
         }

// ---------------------------------------------------------------------------

      } else if ( COLUMN && BACKWARD && LEFT ) {

// ---------------------------------------------------------------------------

         // Let  W =  [ V ]    (M-by-K)
                   // [ I ]    (K-by-K)

         // Form  H C  or  H**T C  where  C = [ B ]  (M-by-N)
                                           // [ A ]  (K-by-N)

         // H = I - W T W**T         or  H**T = I - W T**T W**T

         // A = A -   T (A + V**T B)  or  A = A -   T**T (A + V**T B)
         // B = B - V T (A + V**T B)  or  B = B - V T**T (A + V**T B)

// ---------------------------------------------------------------------------

         MP = min( L+1, M );
         KP = min( K-L+1, K );

         for (J = 1; J <= N; J++) {
            for (I = 1; I <= L; I++) {
               WORK( K-L+I, J ) = B( I, J );
            }
         }

         strmm('L', 'L', 'T', 'N', L, N, ONE, V( 1, KP ), LDV, WORK( KP, 1 ), LDWORK );
         sgemm('T', 'N', L, N, M-L, ONE, V( MP, KP ), LDV, B( MP, 1 ), LDB, ONE, WORK( KP, 1 ), LDWORK );
         sgemm('T', 'N', K-L, N, M, ONE, V, LDV, B, LDB, ZERO, WORK, LDWORK );

         for (J = 1; J <= N; J++) {
            for (I = 1; I <= K; I++) {
               WORK( I, J ) = WORK( I, J ) + A( I, J );
            }
         }

         strmm('L', 'L', TRANS, 'N', K, N, ONE, T, LDT, WORK, LDWORK );

         for (J = 1; J <= N; J++) {
            for (I = 1; I <= K; I++) {
               A( I, J ) = A( I, J ) - WORK( I, J );
            }
         }

         sgemm('N', 'N', M-L, N, K, -ONE, V( MP, 1 ), LDV, WORK, LDWORK, ONE, B( MP, 1 ), LDB );
         sgemm('N', 'N', L, N, K-L, -ONE, V, LDV, WORK, LDWORK, ONE, B,  LDB );
         strmm('L', 'L', 'N', 'N', L, N, ONE, V( 1, KP ), LDV, WORK( KP, 1 ), LDWORK );
         for (J = 1; J <= N; J++) {
            for (I = 1; I <= L; I++) {
               B( I, J ) = B( I, J ) - WORK( K-L+I, J );
            }
         }

// ---------------------------------------------------------------------------

      } else if ( COLUMN && BACKWARD && RIGHT ) {

// ---------------------------------------------------------------------------

         // Let  W =  [ V ]    (N-by-K)
                   // [ I ]    (K-by-K)

         // Form  C H  or  C H**T  where  C = [ B A ] (B is M-by-N, A is M-by-K)

         // H = I - W T W**T          or  H**T = I - W T**T W**T

         // A = A - (A + B V) T       or  A = A - (A + B V) T**T
         // B = B - (A + B V) T V**T  or  B = B - (A + B V) T**T V**T

// ---------------------------------------------------------------------------

         NP = min( L+1, N );
         KP = min( K-L+1, K );

         for (J = 1; J <= L; J++) {
            for (I = 1; I <= M; I++) {
               WORK( I, K-L+J ) = B( I, J );
            }
         }
         strmm('R', 'L', 'N', 'N', M, L, ONE, V( 1, KP ), LDV, WORK( 1, KP ), LDWORK );
         sgemm('N', 'N', M, L, N-L, ONE, B( 1, NP ), LDB, V( NP, KP ), LDV, ONE, WORK( 1, KP ), LDWORK );
         sgemm('N', 'N', M, K-L, N, ONE, B, LDB, V, LDV, ZERO, WORK, LDWORK );

         for (J = 1; J <= K; J++) {
            for (I = 1; I <= M; I++) {
               WORK( I, J ) = WORK( I, J ) + A( I, J );
            }
         }

         strmm('R', 'L', TRANS, 'N', M, K, ONE, T, LDT, WORK, LDWORK );

         for (J = 1; J <= K; J++) {
            for (I = 1; I <= M; I++) {
               A( I, J ) = A( I, J ) - WORK( I, J );
            }
         }

         sgemm('N', 'T', M, N-L, K, -ONE, WORK, LDWORK, V( NP, 1 ), LDV, ONE, B( 1, NP ), LDB );
         sgemm('N', 'T', M, L, K-L, -ONE, WORK, LDWORK, V, LDV, ONE, B, LDB );
         strmm('R', 'L', 'T', 'N', M, L, ONE, V( 1, KP ), LDV, WORK( 1, KP ), LDWORK );
         for (J = 1; J <= L; J++) {
            for (I = 1; I <= M; I++) {
               B( I, J ) = B( I, J ) - WORK( I, K-L+J );
            }
         }

// ---------------------------------------------------------------------------

      } else if ( ROW && FORWARD && LEFT ) {

// ---------------------------------------------------------------------------

         // Let  W =  [ I V ] ( I is K-by-K, V is K-by-M )

         // Form  H C  or  H**T C  where  C = [ A ]  (K-by-N)
                                           // [ B ]  (M-by-N)

         // H = I - W**T T W          or  H**T = I - W**T T**T W

         // A = A -      T (A + V B)  or  A = A -      T**T (A + V B)
         // B = B - V**T T (A + V B)  or  B = B - V**T T**T (A + V B)

// ---------------------------------------------------------------------------

         MP = min( M-L+1, M );
         KP = min( L+1, K );

         for (J = 1; J <= N; J++) {
            for (I = 1; I <= L; I++) {
               WORK( I, J ) = B( M-L+I, J );
            }
         }
         strmm('L', 'L', 'N', 'N', L, N, ONE, V( 1, MP ), LDV, WORK, LDB );
         sgemm('N', 'N', L, N, M-L, ONE, V, LDV,B, LDB, ONE, WORK, LDWORK );
         sgemm('N', 'N', K-L, N, M, ONE, V( KP, 1 ), LDV, B, LDB, ZERO, WORK( KP, 1 ), LDWORK );

         for (J = 1; J <= N; J++) {
            for (I = 1; I <= K; I++) {
               WORK( I, J ) = WORK( I, J ) + A( I, J );
            }
         }

         strmm('L', 'U', TRANS, 'N', K, N, ONE, T, LDT, WORK, LDWORK );

         for (J = 1; J <= N; J++) {
            for (I = 1; I <= K; I++) {
               A( I, J ) = A( I, J ) - WORK( I, J );
            }
         }

         sgemm('T', 'N', M-L, N, K, -ONE, V, LDV, WORK, LDWORK, ONE, B, LDB );
         sgemm('T', 'N', L, N, K-L, -ONE, V( KP, MP ), LDV, WORK( KP, 1 ), LDWORK, ONE, B( MP, 1 ), LDB );
         strmm('L', 'L', 'T', 'N', L, N, ONE, V( 1, MP ), LDV, WORK, LDWORK );
         for (J = 1; J <= N; J++) {
            for (I = 1; I <= L; I++) {
               B( M-L+I, J ) = B( M-L+I, J ) - WORK( I, J );
            }
         }

// ---------------------------------------------------------------------------

      } else if ( ROW && FORWARD && RIGHT ) {

// ---------------------------------------------------------------------------

         // Let  W =  [ I V ] ( I is K-by-K, V is K-by-N )

         // Form  C H  or  C H**T  where  C = [ A B ] (A is M-by-K, B is M-by-N)

         // H = I - W**T T W            or  H**T = I - W**T T**T W

         // A = A - (A + B V**T) T      or  A = A - (A + B V**T) T**T
         // B = B - (A + B V**T) T V    or  B = B - (A + B V**T) T**T V

// ---------------------------------------------------------------------------

         NP = min( N-L+1, N );
         KP = min( L+1, K );

         for (J = 1; J <= L; J++) {
            for (I = 1; I <= M; I++) {
               WORK( I, J ) = B( I, N-L+J );
            }
         }
         strmm('R', 'L', 'T', 'N', M, L, ONE, V( 1, NP ), LDV, WORK, LDWORK );
         sgemm('N', 'T', M, L, N-L, ONE, B, LDB, V, LDV, ONE, WORK, LDWORK );
         sgemm('N', 'T', M, K-L, N, ONE, B, LDB, V( KP, 1 ), LDV, ZERO, WORK( 1, KP ), LDWORK );

         for (J = 1; J <= K; J++) {
            for (I = 1; I <= M; I++) {
               WORK( I, J ) = WORK( I, J ) + A( I, J );
            }
         }

         strmm('R', 'U', TRANS, 'N', M, K, ONE, T, LDT, WORK, LDWORK );

         for (J = 1; J <= K; J++) {
            for (I = 1; I <= M; I++) {
               A( I, J ) = A( I, J ) - WORK( I, J );
            }
         }

         sgemm('N', 'N', M, N-L, K, -ONE, WORK, LDWORK, V, LDV, ONE, B, LDB );
         sgemm('N', 'N', M, L, K-L, -ONE, WORK( 1, KP ), LDWORK, V( KP, NP ), LDV, ONE, B( 1, NP ), LDB );
         strmm('R', 'L', 'N', 'N', M, L, ONE, V( 1, NP ), LDV, WORK, LDWORK );
         for (J = 1; J <= L; J++) {
            for (I = 1; I <= M; I++) {
               B( I, N-L+J ) = B( I, N-L+J ) - WORK( I, J );
            }
         }

// ---------------------------------------------------------------------------

      } else if ( ROW && BACKWARD && LEFT ) {

// ---------------------------------------------------------------------------

         // Let  W =  [ V I ] ( I is K-by-K, V is K-by-M )

         // Form  H C  or  H**T C  where  C = [ B ]  (M-by-N)
                                           // [ A ]  (K-by-N)

         // H = I - W**T T W          or  H**T = I - W**T T**T W

         // A = A -      T (A + V B)  or  A = A -      T**T (A + V B)
         // B = B - V**T T (A + V B)  or  B = B - V**T T**T (A + V B)

// ---------------------------------------------------------------------------

         MP = min( L+1, M );
         KP = min( K-L+1, K );

         for (J = 1; J <= N; J++) {
            for (I = 1; I <= L; I++) {
               WORK( K-L+I, J ) = B( I, J );
            }
         }
         strmm('L', 'U', 'N', 'N', L, N, ONE, V( KP, 1 ), LDV, WORK( KP, 1 ), LDWORK );
         sgemm('N', 'N', L, N, M-L, ONE, V( KP, MP ), LDV, B( MP, 1 ), LDB, ONE, WORK( KP, 1 ), LDWORK );
         sgemm('N', 'N', K-L, N, M, ONE, V, LDV, B, LDB, ZERO, WORK, LDWORK );

         for (J = 1; J <= N; J++) {
            for (I = 1; I <= K; I++) {
               WORK( I, J ) = WORK( I, J ) + A( I, J );
            }
         }

         strmm('L', 'L ', TRANS, 'N', K, N, ONE, T, LDT, WORK, LDWORK );

         for (J = 1; J <= N; J++) {
            for (I = 1; I <= K; I++) {
               A( I, J ) = A( I, J ) - WORK( I, J );
            }
         }

         sgemm('T', 'N', M-L, N, K, -ONE, V( 1, MP ), LDV, WORK, LDWORK, ONE, B( MP, 1 ), LDB );
         sgemm('T', 'N', L, N, K-L, -ONE, V, LDV, WORK, LDWORK, ONE, B, LDB );
         strmm('L', 'U', 'T', 'N', L, N, ONE, V( KP, 1 ), LDV, WORK( KP, 1 ), LDWORK );
         for (J = 1; J <= N; J++) {
            for (I = 1; I <= L; I++) {
               B( I, J ) = B( I, J ) - WORK( K-L+I, J );
            }
         }

// ---------------------------------------------------------------------------

      } else if ( ROW && BACKWARD && RIGHT ) {

// ---------------------------------------------------------------------------

         // Let  W =  [ V I ] ( I is K-by-K, V is K-by-N )

         // Form  C H  or  C H**T  where  C = [ B A ] (A is M-by-K, B is M-by-N)

         // H = I - W**T T W            or  H**T = I - W**T T**T W

         // A = A - (A + B V**T) T      or  A = A - (A + B V**T) T**T
         // B = B - (A + B V**T) T V    or  B = B - (A + B V**T) T**T V

// ---------------------------------------------------------------------------

         NP = min( L+1, N );
         KP = min( K-L+1, K );

         for (J = 1; J <= L; J++) {
            for (I = 1; I <= M; I++) {
               WORK( I, K-L+J ) = B( I, J );
            }
         }
         strmm('R', 'U', 'T', 'N', M, L, ONE, V( KP, 1 ), LDV, WORK( 1, KP ), LDWORK );
         sgemm('N', 'T', M, L, N-L, ONE, B( 1, NP ), LDB, V( KP, NP ), LDV, ONE, WORK( 1, KP ), LDWORK );
         sgemm('N', 'T', M, K-L, N, ONE, B, LDB, V, LDV, ZERO, WORK, LDWORK );

         for (J = 1; J <= K; J++) {
            for (I = 1; I <= M; I++) {
               WORK( I, J ) = WORK( I, J ) + A( I, J );
            }
         }

         strmm('R', 'L', TRANS, 'N', M, K, ONE, T, LDT, WORK, LDWORK );

         for (J = 1; J <= K; J++) {
            for (I = 1; I <= M; I++) {
               A( I, J ) = A( I, J ) - WORK( I, J );
            }
         }

         sgemm('N', 'N', M, N-L, K, -ONE, WORK, LDWORK, V( 1, NP ), LDV, ONE, B( 1, NP ), LDB );
         sgemm('N', 'N', M, L, K-L , -ONE, WORK, LDWORK, V, LDV, ONE, B, LDB );
         strmm('R', 'U', 'N', 'N', M, L, ONE, V( KP, 1 ), LDV, WORK( 1, KP ), LDWORK );
         for (J = 1; J <= L; J++) {
            for (I = 1; I <= M; I++) {
               B( I, J ) = B( I, J ) - WORK( I, K-L+J );
            }
         }

      }

      return;
      }
