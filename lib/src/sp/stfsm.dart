      void stfsm(TRANSR, SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, A, B, LDB ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANSR, DIAG, SIDE, TRANS, UPLO;
      int                LDB, M, N;
      double               ALPHA;
      // ..
      // .. Array Arguments ..
      double               A( 0: * ), B( 0: LDB-1, 0: * );
      // ..

// =====================================================================

      // ..
      // .. Parameters ..
      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      bool               LOWER, LSIDE, MISODD, NISODD, NORMALTRANSR, NOTRANS;
      int                M1, M2, N1, N2, K, INFO, I, J;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, STRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MOD
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      NORMALTRANSR = lsame( TRANSR, 'N' );
      LSIDE = lsame( SIDE, 'L' );
      LOWER = lsame( UPLO, 'L' );
      NOTRANS = lsame( TRANS, 'N' );
      if ( !NORMALTRANSR && !lsame( TRANSR, 'T' ) ) {
         INFO = -1;
      } else if ( !LSIDE && !lsame( SIDE, 'R' ) ) {
         INFO = -2;
      } else if ( !LOWER && !lsame( UPLO, 'U' ) ) {
         INFO = -3;
      } else if ( !NOTRANS && !lsame( TRANS, 'T' ) ) {
         INFO = -4;
      } else if ( !lsame( DIAG, 'N' ) && !lsame( DIAG, 'U' ) ) {
         INFO = -5;
      } else if ( M < 0 ) {
         INFO = -6;
      } else if ( N < 0 ) {
         INFO = -7;
      } else if ( LDB < max( 1, M ) ) {
         INFO = -11;
      }
      if ( INFO != 0 ) {
         xerbla('STFSM ', -INFO );
         return;
      }

      // Quick return when ( (N == 0) || (M == 0) )

      if( ( M == 0 ) || ( N == 0 ) ) return;

      // Quick return when ALPHA == (0D+0)

      if ( ALPHA == ZERO ) {
         for (J = 0; J <= N - 1; J++) { // 20
            for (I = 0; I <= M - 1; I++) { // 10
               B[I, J] = ZERO;
            } // 10
         } // 20
         return;
      }

      if ( LSIDE ) {

         // SIDE = 'L'

         // A is M-by-M.
         // If M is odd, set NISODD = true , and M1 and M2.
         // If M is even, NISODD = false , and M.

         if ( (M % 2) == 0 ) {
            MISODD = false;
            K = M / 2;
         } else {
            MISODD = true;
            if ( LOWER ) {
               M2 = M / 2;
               M1 = M - M2;
            } else {
               M1 = M / 2;
               M2 = M - M1;
            }
         }

         if ( MISODD ) {

            // SIDE = 'L' and N is odd

            if ( NORMALTRANSR ) {

               // SIDE = 'L', N is odd, and TRANSR = 'N'

               if ( LOWER ) {

                  // SIDE  ='L', N is odd, TRANSR = 'N', and UPLO = 'L'

                  if ( NOTRANS ) {

                     // SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'L', and
                     // TRANS = 'N'

                     if ( M == 1 ) {
                        strsm('L', 'L', 'N', DIAG, M1, N, ALPHA, A, M, B, LDB );
                     } else {
                        strsm('L', 'L', 'N', DIAG, M1, N, ALPHA, A( 0 ), M, B, LDB );
                        sgemm('N', 'N', M2, N, M1, -ONE, A( M1 ), M, B, LDB, ALPHA, B( M1, 0 ), LDB );
                        strsm('L', 'U', 'T', DIAG, M2, N, ONE, A( M ), M, B( M1, 0 ), LDB );
                     }

                  } else {

                     // SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'L', and
                     // TRANS = 'T'

                     if ( M == 1 ) {
                        strsm('L', 'L', 'T', DIAG, M1, N, ALPHA, A( 0 ), M, B, LDB );
                     } else {
                        strsm('L', 'U', 'N', DIAG, M2, N, ALPHA, A( M ), M, B( M1, 0 ), LDB );
                        sgemm('T', 'N', M1, N, M2, -ONE, A( M1 ), M, B( M1, 0 ), LDB, ALPHA, B, LDB );
                        strsm('L', 'L', 'T', DIAG, M1, N, ONE, A( 0 ), M, B, LDB );
                     }

                  }

               } else {

                  // SIDE  ='L', N is odd, TRANSR = 'N', and UPLO = 'U'

                  if ( !NOTRANS ) {

                     // SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'U', and
                     // TRANS = 'N'

                     strsm('L', 'L', 'N', DIAG, M1, N, ALPHA, A( M2 ), M, B, LDB );
                     sgemm('T', 'N', M2, N, M1, -ONE, A( 0 ), M, B, LDB, ALPHA, B( M1, 0 ), LDB );
                     strsm('L', 'U', 'T', DIAG, M2, N, ONE, A( M1 ), M, B( M1, 0 ), LDB );

                  } else {

                     // SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'U', and
                     // TRANS = 'T'

                     strsm('L', 'U', 'N', DIAG, M2, N, ALPHA, A( M1 ), M, B( M1, 0 ), LDB );
                     sgemm('N', 'N', M1, N, M2, -ONE, A( 0 ), M, B( M1, 0 ), LDB, ALPHA, B, LDB );
                     strsm('L', 'L', 'T', DIAG, M1, N, ONE, A( M2 ), M, B, LDB );

                  }

               }

            } else {

               // SIDE = 'L', N is odd, and TRANSR = 'T'

               if ( LOWER ) {

                  // SIDE  ='L', N is odd, TRANSR = 'T', and UPLO = 'L'

                  if ( NOTRANS ) {

                     // SIDE  ='L', N is odd, TRANSR = 'T', UPLO = 'L', and
                     // TRANS = 'N'

                     if ( M == 1 ) {
                        strsm('L', 'U', 'T', DIAG, M1, N, ALPHA, A( 0 ), M1, B, LDB );
                     } else {
                        strsm('L', 'U', 'T', DIAG, M1, N, ALPHA, A( 0 ), M1, B, LDB );
                        sgemm('T', 'N', M2, N, M1, -ONE, A( M1*M1 ), M1, B, LDB, ALPHA, B( M1, 0 ), LDB );
                        strsm('L', 'L', 'N', DIAG, M2, N, ONE, A( 1 ), M1, B( M1, 0 ), LDB );
                     }

                  } else {

                     // SIDE  ='L', N is odd, TRANSR = 'T', UPLO = 'L', and
                     // TRANS = 'T'

                     if ( M == 1 ) {
                        strsm('L', 'U', 'N', DIAG, M1, N, ALPHA, A( 0 ), M1, B, LDB );
                     } else {
                        strsm('L', 'L', 'T', DIAG, M2, N, ALPHA, A( 1 ), M1, B( M1, 0 ), LDB );
                        sgemm('N', 'N', M1, N, M2, -ONE, A( M1*M1 ), M1, B( M1, 0 ), LDB, ALPHA, B, LDB );
                        strsm('L', 'U', 'N', DIAG, M1, N, ONE, A( 0 ), M1, B, LDB );
                     }

                  }

               } else {

                  // SIDE  ='L', N is odd, TRANSR = 'T', and UPLO = 'U'

                  if ( !NOTRANS ) {

                     // SIDE  ='L', N is odd, TRANSR = 'T', UPLO = 'U', and
                     // TRANS = 'N'

                     strsm('L', 'U', 'T', DIAG, M1, N, ALPHA, A( M2*M2 ), M2, B, LDB );
                     sgemm('N', 'N', M2, N, M1, -ONE, A( 0 ), M2, B, LDB, ALPHA, B( M1, 0 ), LDB );
                     strsm('L', 'L', 'N', DIAG, M2, N, ONE, A( M1*M2 ), M2, B( M1, 0 ), LDB );

                  } else {

                     // SIDE  ='L', N is odd, TRANSR = 'T', UPLO = 'U', and
                     // TRANS = 'T'

                     strsm('L', 'L', 'T', DIAG, M2, N, ALPHA, A( M1*M2 ), M2, B( M1, 0 ), LDB );
                     sgemm('T', 'N', M1, N, M2, -ONE, A( 0 ), M2, B( M1, 0 ), LDB, ALPHA, B, LDB );
                     strsm('L', 'U', 'N', DIAG, M1, N, ONE, A( M2*M2 ), M2, B, LDB );

                  }

               }

            }

         } else {

            // SIDE = 'L' and N is even

            if ( NORMALTRANSR ) {

               // SIDE = 'L', N is even, and TRANSR = 'N'

               if ( LOWER ) {

                  // SIDE  ='L', N is even, TRANSR = 'N', and UPLO = 'L'

                  if ( NOTRANS ) {

                     // SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'L',
                     // and TRANS = 'N'

                     strsm('L', 'L', 'N', DIAG, K, N, ALPHA, A( 1 ), M+1, B, LDB );
                     sgemm('N', 'N', K, N, K, -ONE, A( K+1 ), M+1, B, LDB, ALPHA, B( K, 0 ), LDB );
                     strsm('L', 'U', 'T', DIAG, K, N, ONE, A( 0 ), M+1, B( K, 0 ), LDB );

                  } else {

                     // SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'L',
                     // and TRANS = 'T'

                     strsm('L', 'U', 'N', DIAG, K, N, ALPHA, A( 0 ), M+1, B( K, 0 ), LDB );
                     sgemm('T', 'N', K, N, K, -ONE, A( K+1 ), M+1, B( K, 0 ), LDB, ALPHA, B, LDB );
                     strsm('L', 'L', 'T', DIAG, K, N, ONE, A( 1 ), M+1, B, LDB );

                  }

               } else {

                  // SIDE  ='L', N is even, TRANSR = 'N', and UPLO = 'U'

                  if ( !NOTRANS ) {

                     // SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'U',
                     // and TRANS = 'N'

                     strsm('L', 'L', 'N', DIAG, K, N, ALPHA, A( K+1 ), M+1, B, LDB );
                     sgemm('T', 'N', K, N, K, -ONE, A( 0 ), M+1, B, LDB, ALPHA, B( K, 0 ), LDB );
                     strsm('L', 'U', 'T', DIAG, K, N, ONE, A( K ), M+1, B( K, 0 ), LDB );

                  } else {

                     // SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'U',
                     // and TRANS = 'T'
                     strsm('L', 'U', 'N', DIAG, K, N, ALPHA, A( K ), M+1, B( K, 0 ), LDB );
                     sgemm('N', 'N', K, N, K, -ONE, A( 0 ), M+1, B( K, 0 ), LDB, ALPHA, B, LDB );
                     strsm('L', 'L', 'T', DIAG, K, N, ONE, A( K+1 ), M+1, B, LDB );

                  }

               }

            } else {

               // SIDE = 'L', N is even, and TRANSR = 'T'

               if ( LOWER ) {

                  // SIDE  ='L', N is even, TRANSR = 'T', and UPLO = 'L'

                  if ( NOTRANS ) {

                     // SIDE  ='L', N is even, TRANSR = 'T', UPLO = 'L',
                     // and TRANS = 'N'

                     strsm('L', 'U', 'T', DIAG, K, N, ALPHA, A( K ), K, B, LDB );
                     sgemm('T', 'N', K, N, K, -ONE, A( K*( K+1 ) ), K, B, LDB, ALPHA, B( K, 0 ), LDB );
                     strsm('L', 'L', 'N', DIAG, K, N, ONE, A( 0 ), K, B( K, 0 ), LDB );

                  } else {

                     // SIDE  ='L', N is even, TRANSR = 'T', UPLO = 'L',
                     // and TRANS = 'T'

                     strsm('L', 'L', 'T', DIAG, K, N, ALPHA, A( 0 ), K, B( K, 0 ), LDB );
                     sgemm('N', 'N', K, N, K, -ONE, A( K*( K+1 ) ), K, B( K, 0 ), LDB, ALPHA, B, LDB );
                     strsm('L', 'U', 'N', DIAG, K, N, ONE, A( K ), K, B, LDB );

                  }

               } else {

                  // SIDE  ='L', N is even, TRANSR = 'T', and UPLO = 'U'

                  if ( !NOTRANS ) {

                     // SIDE  ='L', N is even, TRANSR = 'T', UPLO = 'U',
                     // and TRANS = 'N'

                     strsm('L', 'U', 'T', DIAG, K, N, ALPHA, A( K*( K+1 ) ), K, B, LDB );
                     sgemm('N', 'N', K, N, K, -ONE, A( 0 ), K, B, LDB, ALPHA, B( K, 0 ), LDB );
                     strsm('L', 'L', 'N', DIAG, K, N, ONE, A( K*K ), K, B( K, 0 ), LDB );

                  } else {

                     // SIDE  ='L', N is even, TRANSR = 'T', UPLO = 'U',
                     // and TRANS = 'T'

                     strsm('L', 'L', 'T', DIAG, K, N, ALPHA, A( K*K ), K, B( K, 0 ), LDB );
                     sgemm('T', 'N', K, N, K, -ONE, A( 0 ), K, B( K, 0 ), LDB, ALPHA, B, LDB );
                     strsm('L', 'U', 'N', DIAG, K, N, ONE, A( K*( K+1 ) ), K, B, LDB );

                  }

               }

            }

         }

      } else {

         // SIDE = 'R'

         // A is N-by-N.
         // If N is odd, set NISODD = true , and N1 and N2.
         // If N is even, NISODD = false , and K.

         if ( (N % 2) == 0 ) {
            NISODD = false;
            K = N / 2;
         } else {
            NISODD = true;
            if ( LOWER ) {
               N2 = N / 2;
               N1 = N - N2;
            } else {
               N1 = N / 2;
               N2 = N - N1;
            }
         }

         if ( NISODD ) {

            // SIDE = 'R' and N is odd

            if ( NORMALTRANSR ) {

               // SIDE = 'R', N is odd, and TRANSR = 'N'

               if ( LOWER ) {

                  // SIDE  ='R', N is odd, TRANSR = 'N', and UPLO = 'L'

                  if ( NOTRANS ) {

                     // SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'L', and
                     // TRANS = 'N'

                     strsm('R', 'U', 'T', DIAG, M, N2, ALPHA, A( N ), N, B( 0, N1 ), LDB );
                     sgemm('N', 'N', M, N1, N2, -ONE, B( 0, N1 ), LDB, A( N1 ), N, ALPHA, B( 0, 0 ), LDB );
                     strsm('R', 'L', 'N', DIAG, M, N1, ONE, A( 0 ), N, B( 0, 0 ), LDB );

                  } else {

                     // SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'L', and
                     // TRANS = 'T'

                     strsm('R', 'L', 'T', DIAG, M, N1, ALPHA, A( 0 ), N, B( 0, 0 ), LDB );
                     sgemm('N', 'T', M, N2, N1, -ONE, B( 0, 0 ), LDB, A( N1 ), N, ALPHA, B( 0, N1 ), LDB );
                     strsm('R', 'U', 'N', DIAG, M, N2, ONE, A( N ), N, B( 0, N1 ), LDB );

                  }

               } else {

                  // SIDE  ='R', N is odd, TRANSR = 'N', and UPLO = 'U'

                  if ( NOTRANS ) {

                     // SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'U', and
                     // TRANS = 'N'

                     strsm('R', 'L', 'T', DIAG, M, N1, ALPHA, A( N2 ), N, B( 0, 0 ), LDB );
                     sgemm('N', 'N', M, N2, N1, -ONE, B( 0, 0 ), LDB, A( 0 ), N, ALPHA, B( 0, N1 ), LDB );
                     strsm('R', 'U', 'N', DIAG, M, N2, ONE, A( N1 ), N, B( 0, N1 ), LDB );

                  } else {

                     // SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'U', and
                     // TRANS = 'T'

                     strsm('R', 'U', 'T', DIAG, M, N2, ALPHA, A( N1 ), N, B( 0, N1 ), LDB );
                     sgemm('N', 'T', M, N1, N2, -ONE, B( 0, N1 ), LDB, A( 0 ), N, ALPHA, B( 0, 0 ), LDB );
                     strsm('R', 'L', 'N', DIAG, M, N1, ONE, A( N2 ), N, B( 0, 0 ), LDB );

                  }

               }

            } else {

               // SIDE = 'R', N is odd, and TRANSR = 'T'

               if ( LOWER ) {

                  // SIDE  ='R', N is odd, TRANSR = 'T', and UPLO = 'L'

                  if ( NOTRANS ) {

                     // SIDE  ='R', N is odd, TRANSR = 'T', UPLO = 'L', and
                     // TRANS = 'N'

                     strsm('R', 'L', 'N', DIAG, M, N2, ALPHA, A( 1 ), N1, B( 0, N1 ), LDB );
                     sgemm('N', 'T', M, N1, N2, -ONE, B( 0, N1 ), LDB, A( N1*N1 ), N1, ALPHA, B( 0, 0 ), LDB );
                     strsm('R', 'U', 'T', DIAG, M, N1, ONE, A( 0 ), N1, B( 0, 0 ), LDB );

                  } else {

                     // SIDE  ='R', N is odd, TRANSR = 'T', UPLO = 'L', and
                     // TRANS = 'T'

                     strsm('R', 'U', 'N', DIAG, M, N1, ALPHA, A( 0 ), N1, B( 0, 0 ), LDB );
                     sgemm('N', 'N', M, N2, N1, -ONE, B( 0, 0 ), LDB, A( N1*N1 ), N1, ALPHA, B( 0, N1 ), LDB );
                     strsm('R', 'L', 'T', DIAG, M, N2, ONE, A( 1 ), N1, B( 0, N1 ), LDB );

                  }

               } else {

                  // SIDE  ='R', N is odd, TRANSR = 'T', and UPLO = 'U'

                  if ( NOTRANS ) {

                     // SIDE  ='R', N is odd, TRANSR = 'T', UPLO = 'U', and
                     // TRANS = 'N'

                     strsm('R', 'U', 'N', DIAG, M, N1, ALPHA, A( N2*N2 ), N2, B( 0, 0 ), LDB );
                     sgemm('N', 'T', M, N2, N1, -ONE, B( 0, 0 ), LDB, A( 0 ), N2, ALPHA, B( 0, N1 ), LDB );
                     strsm('R', 'L', 'T', DIAG, M, N2, ONE, A( N1*N2 ), N2, B( 0, N1 ), LDB );

                  } else {

                     // SIDE  ='R', N is odd, TRANSR = 'T', UPLO = 'U', and
                     // TRANS = 'T'

                     strsm('R', 'L', 'N', DIAG, M, N2, ALPHA, A( N1*N2 ), N2, B( 0, N1 ), LDB );
                     sgemm('N', 'N', M, N1, N2, -ONE, B( 0, N1 ), LDB, A( 0 ), N2, ALPHA, B( 0, 0 ), LDB );
                     strsm('R', 'U', 'T', DIAG, M, N1, ONE, A( N2*N2 ), N2, B( 0, 0 ), LDB );

                  }

               }

            }

         } else {

            // SIDE = 'R' and N is even

            if ( NORMALTRANSR ) {

               // SIDE = 'R', N is even, and TRANSR = 'N'

               if ( LOWER ) {

                  // SIDE  ='R', N is even, TRANSR = 'N', and UPLO = 'L'

                  if ( NOTRANS ) {

                     // SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'L',
                     // and TRANS = 'N'

                     strsm('R', 'U', 'T', DIAG, M, K, ALPHA, A( 0 ), N+1, B( 0, K ), LDB );
                     sgemm('N', 'N', M, K, K, -ONE, B( 0, K ), LDB, A( K+1 ), N+1, ALPHA, B( 0, 0 ), LDB );
                     strsm('R', 'L', 'N', DIAG, M, K, ONE, A( 1 ), N+1, B( 0, 0 ), LDB );

                  } else {

                     // SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'L',
                     // and TRANS = 'T'

                     strsm('R', 'L', 'T', DIAG, M, K, ALPHA, A( 1 ), N+1, B( 0, 0 ), LDB );
                     sgemm('N', 'T', M, K, K, -ONE, B( 0, 0 ), LDB, A( K+1 ), N+1, ALPHA, B( 0, K ), LDB );
                     strsm('R', 'U', 'N', DIAG, M, K, ONE, A( 0 ), N+1, B( 0, K ), LDB );

                  }

               } else {

                  // SIDE  ='R', N is even, TRANSR = 'N', and UPLO = 'U'

                  if ( NOTRANS ) {

                     // SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'U',
                     // and TRANS = 'N'

                     strsm('R', 'L', 'T', DIAG, M, K, ALPHA, A( K+1 ), N+1, B( 0, 0 ), LDB );
                     sgemm('N', 'N', M, K, K, -ONE, B( 0, 0 ), LDB, A( 0 ), N+1, ALPHA, B( 0, K ), LDB );
                     strsm('R', 'U', 'N', DIAG, M, K, ONE, A( K ), N+1, B( 0, K ), LDB );

                  } else {

                     // SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'U',
                     // and TRANS = 'T'

                     strsm('R', 'U', 'T', DIAG, M, K, ALPHA, A( K ), N+1, B( 0, K ), LDB );
                     sgemm('N', 'T', M, K, K, -ONE, B( 0, K ), LDB, A( 0 ), N+1, ALPHA, B( 0, 0 ), LDB );
                     strsm('R', 'L', 'N', DIAG, M, K, ONE, A( K+1 ), N+1, B( 0, 0 ), LDB );

                  }

               }

            } else {

               // SIDE = 'R', N is even, and TRANSR = 'T'

               if ( LOWER ) {

                  // SIDE  ='R', N is even, TRANSR = 'T', and UPLO = 'L'

                  if ( NOTRANS ) {

                     // SIDE  ='R', N is even, TRANSR = 'T', UPLO = 'L',
                     // and TRANS = 'N'

                     strsm('R', 'L', 'N', DIAG, M, K, ALPHA, A( 0 ), K, B( 0, K ), LDB );
                     sgemm('N', 'T', M, K, K, -ONE, B( 0, K ), LDB, A( ( K+1 )*K ), K, ALPHA, B( 0, 0 ), LDB );
                     strsm('R', 'U', 'T', DIAG, M, K, ONE, A( K ), K, B( 0, 0 ), LDB );

                  } else {

                     // SIDE  ='R', N is even, TRANSR = 'T', UPLO = 'L',
                     // and TRANS = 'T'

                     strsm('R', 'U', 'N', DIAG, M, K, ALPHA, A( K ), K, B( 0, 0 ), LDB );
                     sgemm('N', 'N', M, K, K, -ONE, B( 0, 0 ), LDB, A( ( K+1 )*K ), K, ALPHA, B( 0, K ), LDB );
                     strsm('R', 'L', 'T', DIAG, M, K, ONE, A( 0 ), K, B( 0, K ), LDB );

                  }

               } else {

                  // SIDE  ='R', N is even, TRANSR = 'T', and UPLO = 'U'

                  if ( NOTRANS ) {

                     // SIDE  ='R', N is even, TRANSR = 'T', UPLO = 'U',
                     // and TRANS = 'N'

                     strsm('R', 'U', 'N', DIAG, M, K, ALPHA, A( ( K+1 )*K ), K, B( 0, 0 ), LDB );
                     sgemm('N', 'T', M, K, K, -ONE, B( 0, 0 ), LDB, A( 0 ), K, ALPHA, B( 0, K ), LDB );
                     strsm('R', 'L', 'T', DIAG, M, K, ONE, A( K*K ), K, B( 0, K ), LDB );

                  } else {

                     // SIDE  ='R', N is even, TRANSR = 'T', UPLO = 'U',
                     // and TRANS = 'T'

                     strsm('R', 'L', 'N', DIAG, M, K, ALPHA, A( K*K ), K, B( 0, K ), LDB );
                     sgemm('N', 'N', M, K, K, -ONE, B( 0, K ), LDB, A( 0 ), K, ALPHA, B( 0, 0 ), LDB );
                     strsm('R', 'U', 'T', DIAG, M, K, ONE, A( ( K+1 )*K ), K, B( 0, 0 ), LDB );

                  }

               }

            }

         }
      }

      return;
      }