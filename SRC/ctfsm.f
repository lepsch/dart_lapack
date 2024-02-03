      SUBROUTINE CTFSM( TRANSR, SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, A, B, LDB );

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANSR, DIAG, SIDE, TRANS, UPLO;
      int                LDB, M, N;
      COMPLEX            ALPHA;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( 0: * ), B( 0: LDB-1, 0: * );
      // ..

*  =====================================================================
      // ..
      // .. Parameters ..
      COMPLEX            CONE, CZERO;
      const              CONE = ( 1.0, 0.0 ), CZERO = ( 0.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LOWER, LSIDE, MISODD, NISODD, NORMALTRANSR, NOTRANS;
      int                M1, M2, N1, N2, K, INFO, I, J;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, CGEMM, CTRSM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MOD
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      NORMALTRANSR = LSAME( TRANSR, 'N' );
      LSIDE = LSAME( SIDE, 'L' );
      LOWER = LSAME( UPLO, 'L' );
      NOTRANS = LSAME( TRANS, 'N' );
      if ( !NORMALTRANSR && !LSAME( TRANSR, 'C' ) ) {
         INFO = -1;
      } else if ( !LSIDE && !LSAME( SIDE, 'R' ) ) {
         INFO = -2;
      } else if ( !LOWER && !LSAME( UPLO, 'U' ) ) {
         INFO = -3;
      } else if ( !NOTRANS && !LSAME( TRANS, 'C' ) ) {
         INFO = -4;
      } else if ( !LSAME( DIAG, 'N' ) && !LSAME( DIAG, 'U' ) ) {
         INFO = -5;
      } else if ( M < 0 ) {
         INFO = -6;
      } else if ( N < 0 ) {
         INFO = -7;
      } else if ( LDB < MAX( 1, M ) ) {
         INFO = -11;
      }
      if ( INFO != 0 ) {
         xerbla('CTFSM ', -INFO );
         RETURN;
      }

      // Quick return when ( (N == 0) || (M == 0) )

      IF( ( M == 0 ) || ( N == 0 ) ) RETURN;

      // Quick return when ALPHA == (0E+0,0E+0)

      if ( ALPHA == CZERO ) {
         for (J = 0; J <= N - 1; J++) { // 20
            for (I = 0; I <= M - 1; I++) { // 10
               B( I, J ) = CZERO;
            } // 10
         } // 20
         RETURN;
      }

      if ( LSIDE ) {

         // SIDE = 'L'

         // A is M-by-M.
         // If M is odd, set NISODD = true , and M1 and M2.
         // If M is even, NISODD = false , and M.

         if ( MOD( M, 2 ) == 0 ) {
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
                        ctrsm('L', 'L', 'N', DIAG, M1, N, ALPHA, A, M, B, LDB );
                     } else {
                        ctrsm('L', 'L', 'N', DIAG, M1, N, ALPHA, A( 0 ), M, B, LDB );
                        cgemm('N', 'N', M2, N, M1, -CONE, A( M1 ), M, B, LDB, ALPHA, B( M1, 0 ), LDB );
                        ctrsm('L', 'U', 'C', DIAG, M2, N, CONE, A( M ), M, B( M1, 0 ), LDB );
                     }

                  } else {

                     // SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'L', and
                     // TRANS = 'C'

                     if ( M == 1 ) {
                        ctrsm('L', 'L', 'C', DIAG, M1, N, ALPHA, A( 0 ), M, B, LDB );
                     } else {
                        ctrsm('L', 'U', 'N', DIAG, M2, N, ALPHA, A( M ), M, B( M1, 0 ), LDB );
                        cgemm('C', 'N', M1, N, M2, -CONE, A( M1 ), M, B( M1, 0 ), LDB, ALPHA, B, LDB );
                        ctrsm('L', 'L', 'C', DIAG, M1, N, CONE, A( 0 ), M, B, LDB );
                     }

                  }

               } else {

                  // SIDE  ='L', N is odd, TRANSR = 'N', and UPLO = 'U'

                  if ( !NOTRANS ) {

                     // SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'U', and
                     // TRANS = 'N'

                     ctrsm('L', 'L', 'N', DIAG, M1, N, ALPHA, A( M2 ), M, B, LDB );
                     cgemm('C', 'N', M2, N, M1, -CONE, A( 0 ), M, B, LDB, ALPHA, B( M1, 0 ), LDB );
                     ctrsm('L', 'U', 'C', DIAG, M2, N, CONE, A( M1 ), M, B( M1, 0 ), LDB );

                  } else {

                     // SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'U', and
                     // TRANS = 'C'

                     ctrsm('L', 'U', 'N', DIAG, M2, N, ALPHA, A( M1 ), M, B( M1, 0 ), LDB );
                     cgemm('N', 'N', M1, N, M2, -CONE, A( 0 ), M, B( M1, 0 ), LDB, ALPHA, B, LDB );
                     ctrsm('L', 'L', 'C', DIAG, M1, N, CONE, A( M2 ), M, B, LDB );

                  }

               }

            } else {

               // SIDE = 'L', N is odd, and TRANSR = 'C'

               if ( LOWER ) {

                  // SIDE  ='L', N is odd, TRANSR = 'C', and UPLO = 'L'

                  if ( NOTRANS ) {

                     // SIDE  ='L', N is odd, TRANSR = 'C', UPLO = 'L', and
                     // TRANS = 'N'

                     if ( M == 1 ) {
                        ctrsm('L', 'U', 'C', DIAG, M1, N, ALPHA, A( 0 ), M1, B, LDB );
                     } else {
                        ctrsm('L', 'U', 'C', DIAG, M1, N, ALPHA, A( 0 ), M1, B, LDB );
                        cgemm('C', 'N', M2, N, M1, -CONE, A( M1*M1 ), M1, B, LDB, ALPHA, B( M1, 0 ), LDB );
                        ctrsm('L', 'L', 'N', DIAG, M2, N, CONE, A( 1 ), M1, B( M1, 0 ), LDB );
                     }

                  } else {

                     // SIDE  ='L', N is odd, TRANSR = 'C', UPLO = 'L', and
                     // TRANS = 'C'

                     if ( M == 1 ) {
                        ctrsm('L', 'U', 'N', DIAG, M1, N, ALPHA, A( 0 ), M1, B, LDB );
                     } else {
                        ctrsm('L', 'L', 'C', DIAG, M2, N, ALPHA, A( 1 ), M1, B( M1, 0 ), LDB );
                        cgemm('N', 'N', M1, N, M2, -CONE, A( M1*M1 ), M1, B( M1, 0 ), LDB, ALPHA, B, LDB );
                        ctrsm('L', 'U', 'N', DIAG, M1, N, CONE, A( 0 ), M1, B, LDB );
                     }

                  }

               } else {

                  // SIDE  ='L', N is odd, TRANSR = 'C', and UPLO = 'U'

                  if ( !NOTRANS ) {

                     // SIDE  ='L', N is odd, TRANSR = 'C', UPLO = 'U', and
                     // TRANS = 'N'

                     ctrsm('L', 'U', 'C', DIAG, M1, N, ALPHA, A( M2*M2 ), M2, B, LDB );
                     cgemm('N', 'N', M2, N, M1, -CONE, A( 0 ), M2, B, LDB, ALPHA, B( M1, 0 ), LDB );
                     ctrsm('L', 'L', 'N', DIAG, M2, N, CONE, A( M1*M2 ), M2, B( M1, 0 ), LDB );

                  } else {

                     // SIDE  ='L', N is odd, TRANSR = 'C', UPLO = 'U', and
                     // TRANS = 'C'

                     ctrsm('L', 'L', 'C', DIAG, M2, N, ALPHA, A( M1*M2 ), M2, B( M1, 0 ), LDB );
                     cgemm('C', 'N', M1, N, M2, -CONE, A( 0 ), M2, B( M1, 0 ), LDB, ALPHA, B, LDB );
                     ctrsm('L', 'U', 'N', DIAG, M1, N, CONE, A( M2*M2 ), M2, B, LDB );

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

                     ctrsm('L', 'L', 'N', DIAG, K, N, ALPHA, A( 1 ), M+1, B, LDB );
                     cgemm('N', 'N', K, N, K, -CONE, A( K+1 ), M+1, B, LDB, ALPHA, B( K, 0 ), LDB );
                     ctrsm('L', 'U', 'C', DIAG, K, N, CONE, A( 0 ), M+1, B( K, 0 ), LDB );

                  } else {

                     // SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'L',
                     // and TRANS = 'C'

                     ctrsm('L', 'U', 'N', DIAG, K, N, ALPHA, A( 0 ), M+1, B( K, 0 ), LDB );
                     cgemm('C', 'N', K, N, K, -CONE, A( K+1 ), M+1, B( K, 0 ), LDB, ALPHA, B, LDB );
                     ctrsm('L', 'L', 'C', DIAG, K, N, CONE, A( 1 ), M+1, B, LDB );

                  }

               } else {

                  // SIDE  ='L', N is even, TRANSR = 'N', and UPLO = 'U'

                  if ( !NOTRANS ) {

                     // SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'U',
                     // and TRANS = 'N'

                     ctrsm('L', 'L', 'N', DIAG, K, N, ALPHA, A( K+1 ), M+1, B, LDB );
                     cgemm('C', 'N', K, N, K, -CONE, A( 0 ), M+1, B, LDB, ALPHA, B( K, 0 ), LDB );
                     ctrsm('L', 'U', 'C', DIAG, K, N, CONE, A( K ), M+1, B( K, 0 ), LDB );

                  } else {

                     // SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'U',
                     // and TRANS = 'C'
                     ctrsm('L', 'U', 'N', DIAG, K, N, ALPHA, A( K ), M+1, B( K, 0 ), LDB );
                     cgemm('N', 'N', K, N, K, -CONE, A( 0 ), M+1, B( K, 0 ), LDB, ALPHA, B, LDB );
                     ctrsm('L', 'L', 'C', DIAG, K, N, CONE, A( K+1 ), M+1, B, LDB );

                  }

               }

            } else {

               // SIDE = 'L', N is even, and TRANSR = 'C'

               if ( LOWER ) {

                  // SIDE  ='L', N is even, TRANSR = 'C', and UPLO = 'L'

                  if ( NOTRANS ) {

                     // SIDE  ='L', N is even, TRANSR = 'C', UPLO = 'L',
                     // and TRANS = 'N'

                     ctrsm('L', 'U', 'C', DIAG, K, N, ALPHA, A( K ), K, B, LDB );
                     cgemm('C', 'N', K, N, K, -CONE, A( K*( K+1 ) ), K, B, LDB, ALPHA, B( K, 0 ), LDB );
                     ctrsm('L', 'L', 'N', DIAG, K, N, CONE, A( 0 ), K, B( K, 0 ), LDB );

                  } else {

                     // SIDE  ='L', N is even, TRANSR = 'C', UPLO = 'L',
                     // and TRANS = 'C'

                     ctrsm('L', 'L', 'C', DIAG, K, N, ALPHA, A( 0 ), K, B( K, 0 ), LDB );
                     cgemm('N', 'N', K, N, K, -CONE, A( K*( K+1 ) ), K, B( K, 0 ), LDB, ALPHA, B, LDB );
                     ctrsm('L', 'U', 'N', DIAG, K, N, CONE, A( K ), K, B, LDB );

                  }

               } else {

                  // SIDE  ='L', N is even, TRANSR = 'C', and UPLO = 'U'

                  if ( !NOTRANS ) {

                     // SIDE  ='L', N is even, TRANSR = 'C', UPLO = 'U',
                     // and TRANS = 'N'

                     ctrsm('L', 'U', 'C', DIAG, K, N, ALPHA, A( K*( K+1 ) ), K, B, LDB );
                     cgemm('N', 'N', K, N, K, -CONE, A( 0 ), K, B, LDB, ALPHA, B( K, 0 ), LDB );
                     ctrsm('L', 'L', 'N', DIAG, K, N, CONE, A( K*K ), K, B( K, 0 ), LDB );

                  } else {

                     // SIDE  ='L', N is even, TRANSR = 'C', UPLO = 'U',
                     // and TRANS = 'C'

                     ctrsm('L', 'L', 'C', DIAG, K, N, ALPHA, A( K*K ), K, B( K, 0 ), LDB );
                     cgemm('C', 'N', K, N, K, -CONE, A( 0 ), K, B( K, 0 ), LDB, ALPHA, B, LDB );
                     ctrsm('L', 'U', 'N', DIAG, K, N, CONE, A( K*( K+1 ) ), K, B, LDB );

                  }

               }

            }

         }

      } else {

         // SIDE = 'R'

         // A is N-by-N.
         // If N is odd, set NISODD = true , and N1 and N2.
         // If N is even, NISODD = false , and K.

         if ( MOD( N, 2 ) == 0 ) {
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

                     ctrsm('R', 'U', 'C', DIAG, M, N2, ALPHA, A( N ), N, B( 0, N1 ), LDB );
                     cgemm('N', 'N', M, N1, N2, -CONE, B( 0, N1 ), LDB, A( N1 ), N, ALPHA, B( 0, 0 ), LDB );
                     ctrsm('R', 'L', 'N', DIAG, M, N1, CONE, A( 0 ), N, B( 0, 0 ), LDB );

                  } else {

                     // SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'L', and
                     // TRANS = 'C'

                     ctrsm('R', 'L', 'C', DIAG, M, N1, ALPHA, A( 0 ), N, B( 0, 0 ), LDB );
                     cgemm('N', 'C', M, N2, N1, -CONE, B( 0, 0 ), LDB, A( N1 ), N, ALPHA, B( 0, N1 ), LDB );
                     ctrsm('R', 'U', 'N', DIAG, M, N2, CONE, A( N ), N, B( 0, N1 ), LDB );

                  }

               } else {

                  // SIDE  ='R', N is odd, TRANSR = 'N', and UPLO = 'U'

                  if ( NOTRANS ) {

                     // SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'U', and
                     // TRANS = 'N'

                     ctrsm('R', 'L', 'C', DIAG, M, N1, ALPHA, A( N2 ), N, B( 0, 0 ), LDB );
                     cgemm('N', 'N', M, N2, N1, -CONE, B( 0, 0 ), LDB, A( 0 ), N, ALPHA, B( 0, N1 ), LDB );
                     ctrsm('R', 'U', 'N', DIAG, M, N2, CONE, A( N1 ), N, B( 0, N1 ), LDB );

                  } else {

                     // SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'U', and
                     // TRANS = 'C'

                     ctrsm('R', 'U', 'C', DIAG, M, N2, ALPHA, A( N1 ), N, B( 0, N1 ), LDB );
                     cgemm('N', 'C', M, N1, N2, -CONE, B( 0, N1 ), LDB, A( 0 ), N, ALPHA, B( 0, 0 ), LDB );
                     ctrsm('R', 'L', 'N', DIAG, M, N1, CONE, A( N2 ), N, B( 0, 0 ), LDB );

                  }

               }

            } else {

               // SIDE = 'R', N is odd, and TRANSR = 'C'

               if ( LOWER ) {

                  // SIDE  ='R', N is odd, TRANSR = 'C', and UPLO = 'L'

                  if ( NOTRANS ) {

                     // SIDE  ='R', N is odd, TRANSR = 'C', UPLO = 'L', and
                     // TRANS = 'N'

                     ctrsm('R', 'L', 'N', DIAG, M, N2, ALPHA, A( 1 ), N1, B( 0, N1 ), LDB );
                     cgemm('N', 'C', M, N1, N2, -CONE, B( 0, N1 ), LDB, A( N1*N1 ), N1, ALPHA, B( 0, 0 ), LDB );
                     ctrsm('R', 'U', 'C', DIAG, M, N1, CONE, A( 0 ), N1, B( 0, 0 ), LDB );

                  } else {

                     // SIDE  ='R', N is odd, TRANSR = 'C', UPLO = 'L', and
                     // TRANS = 'C'

                     ctrsm('R', 'U', 'N', DIAG, M, N1, ALPHA, A( 0 ), N1, B( 0, 0 ), LDB );
                     cgemm('N', 'N', M, N2, N1, -CONE, B( 0, 0 ), LDB, A( N1*N1 ), N1, ALPHA, B( 0, N1 ), LDB );
                     ctrsm('R', 'L', 'C', DIAG, M, N2, CONE, A( 1 ), N1, B( 0, N1 ), LDB );

                  }

               } else {

                  // SIDE  ='R', N is odd, TRANSR = 'C', and UPLO = 'U'

                  if ( NOTRANS ) {

                     // SIDE  ='R', N is odd, TRANSR = 'C', UPLO = 'U', and
                     // TRANS = 'N'

                     ctrsm('R', 'U', 'N', DIAG, M, N1, ALPHA, A( N2*N2 ), N2, B( 0, 0 ), LDB );
                     cgemm('N', 'C', M, N2, N1, -CONE, B( 0, 0 ), LDB, A( 0 ), N2, ALPHA, B( 0, N1 ), LDB );
                     ctrsm('R', 'L', 'C', DIAG, M, N2, CONE, A( N1*N2 ), N2, B( 0, N1 ), LDB );

                  } else {

                     // SIDE  ='R', N is odd, TRANSR = 'C', UPLO = 'U', and
                     // TRANS = 'C'

                     ctrsm('R', 'L', 'N', DIAG, M, N2, ALPHA, A( N1*N2 ), N2, B( 0, N1 ), LDB );
                     cgemm('N', 'N', M, N1, N2, -CONE, B( 0, N1 ), LDB, A( 0 ), N2, ALPHA, B( 0, 0 ), LDB );
                     ctrsm('R', 'U', 'C', DIAG, M, N1, CONE, A( N2*N2 ), N2, B( 0, 0 ), LDB );

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

                     ctrsm('R', 'U', 'C', DIAG, M, K, ALPHA, A( 0 ), N+1, B( 0, K ), LDB );
                     cgemm('N', 'N', M, K, K, -CONE, B( 0, K ), LDB, A( K+1 ), N+1, ALPHA, B( 0, 0 ), LDB );
                     ctrsm('R', 'L', 'N', DIAG, M, K, CONE, A( 1 ), N+1, B( 0, 0 ), LDB );

                  } else {

                     // SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'L',
                     // and TRANS = 'C'

                     ctrsm('R', 'L', 'C', DIAG, M, K, ALPHA, A( 1 ), N+1, B( 0, 0 ), LDB );
                     cgemm('N', 'C', M, K, K, -CONE, B( 0, 0 ), LDB, A( K+1 ), N+1, ALPHA, B( 0, K ), LDB );
                     ctrsm('R', 'U', 'N', DIAG, M, K, CONE, A( 0 ), N+1, B( 0, K ), LDB );

                  }

               } else {

                  // SIDE  ='R', N is even, TRANSR = 'N', and UPLO = 'U'

                  if ( NOTRANS ) {

                     // SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'U',
                     // and TRANS = 'N'

                     ctrsm('R', 'L', 'C', DIAG, M, K, ALPHA, A( K+1 ), N+1, B( 0, 0 ), LDB );
                     cgemm('N', 'N', M, K, K, -CONE, B( 0, 0 ), LDB, A( 0 ), N+1, ALPHA, B( 0, K ), LDB );
                     ctrsm('R', 'U', 'N', DIAG, M, K, CONE, A( K ), N+1, B( 0, K ), LDB );

                  } else {

                     // SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'U',
                     // and TRANS = 'C'

                     ctrsm('R', 'U', 'C', DIAG, M, K, ALPHA, A( K ), N+1, B( 0, K ), LDB );
                     cgemm('N', 'C', M, K, K, -CONE, B( 0, K ), LDB, A( 0 ), N+1, ALPHA, B( 0, 0 ), LDB );
                     ctrsm('R', 'L', 'N', DIAG, M, K, CONE, A( K+1 ), N+1, B( 0, 0 ), LDB );

                  }

               }

            } else {

               // SIDE = 'R', N is even, and TRANSR = 'C'

               if ( LOWER ) {

                  // SIDE  ='R', N is even, TRANSR = 'C', and UPLO = 'L'

                  if ( NOTRANS ) {

                     // SIDE  ='R', N is even, TRANSR = 'C', UPLO = 'L',
                     // and TRANS = 'N'

                     ctrsm('R', 'L', 'N', DIAG, M, K, ALPHA, A( 0 ), K, B( 0, K ), LDB );
                     cgemm('N', 'C', M, K, K, -CONE, B( 0, K ), LDB, A( ( K+1 )*K ), K, ALPHA, B( 0, 0 ), LDB );
                     ctrsm('R', 'U', 'C', DIAG, M, K, CONE, A( K ), K, B( 0, 0 ), LDB );

                  } else {

                     // SIDE  ='R', N is even, TRANSR = 'C', UPLO = 'L',
                     // and TRANS = 'C'

                     ctrsm('R', 'U', 'N', DIAG, M, K, ALPHA, A( K ), K, B( 0, 0 ), LDB );
                     cgemm('N', 'N', M, K, K, -CONE, B( 0, 0 ), LDB, A( ( K+1 )*K ), K, ALPHA, B( 0, K ), LDB );
                     ctrsm('R', 'L', 'C', DIAG, M, K, CONE, A( 0 ), K, B( 0, K ), LDB );

                  }

               } else {

                  // SIDE  ='R', N is even, TRANSR = 'C', and UPLO = 'U'

                  if ( NOTRANS ) {

                     // SIDE  ='R', N is even, TRANSR = 'C', UPLO = 'U',
                     // and TRANS = 'N'

                     ctrsm('R', 'U', 'N', DIAG, M, K, ALPHA, A( ( K+1 )*K ), K, B( 0, 0 ), LDB );
                     cgemm('N', 'C', M, K, K, -CONE, B( 0, 0 ), LDB, A( 0 ), K, ALPHA, B( 0, K ), LDB );
                     ctrsm('R', 'L', 'C', DIAG, M, K, CONE, A( K*K ), K, B( 0, K ), LDB );

                  } else {

                     // SIDE  ='R', N is even, TRANSR = 'C', UPLO = 'U',
                     // and TRANS = 'C'

                     ctrsm('R', 'L', 'N', DIAG, M, K, ALPHA, A( K*K ), K, B( 0, K ), LDB );
                     cgemm('N', 'N', M, K, K, -CONE, B( 0, K ), LDB, A( 0 ), K, ALPHA, B( 0, 0 ), LDB );
                     ctrsm('R', 'U', 'C', DIAG, M, K, CONE, A( ( K+1 )*K ), K, B( 0, 0 ), LDB );

                  }

               }

            }

         }
      }

      RETURN;

      // End of CTFSM

      }
