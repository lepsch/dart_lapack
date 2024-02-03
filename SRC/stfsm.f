      SUBROUTINE STFSM( TRANSR, SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, A, B, LDB )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANSR, DIAG, SIDE, TRANS, UPLO;
      int                LDB, M, N;
      REAL               ALPHA
      // ..
      // .. Array Arguments ..
      REAL               A( 0: * ), B( 0: LDB-1, 0: * )
      // ..

*  =====================================================================

      // ..
      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
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
      // EXTERNAL SGEMM, STRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MOD
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      NORMALTRANSR = LSAME( TRANSR, 'N' )
      LSIDE = LSAME( SIDE, 'L' )
      LOWER = LSAME( UPLO, 'L' )
      NOTRANS = LSAME( TRANS, 'N' )
      if ( .NOT.NORMALTRANSR .AND. .NOT.LSAME( TRANSR, 'T' ) ) {
         INFO = -1
      } else if ( .NOT.LSIDE .AND. .NOT.LSAME( SIDE, 'R' ) ) {
         INFO = -2
      } else if ( .NOT.LOWER .AND. .NOT.LSAME( UPLO, 'U' ) ) {
         INFO = -3
      } else if ( .NOT.NOTRANS .AND. .NOT.LSAME( TRANS, 'T' ) ) {
         INFO = -4
      } else if ( .NOT.LSAME( DIAG, 'N' ) .AND. .NOT.LSAME( DIAG, 'U' ) ) {
         INFO = -5
      } else if ( M.LT.0 ) {
         INFO = -6
      } else if ( N.LT.0 ) {
         INFO = -7
      } else if ( LDB.LT.MAX( 1, M ) ) {
         INFO = -11
      }
      if ( INFO.NE.0 ) {
         xerbla('STFSM ', -INFO );
         RETURN
      }

      // Quick return when ( (N.EQ.0).OR.(M.EQ.0) )

      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) ) RETURN

      // Quick return when ALPHA.EQ.(0D+0)

      if ( ALPHA.EQ.ZERO ) {
         DO 20 J = 0, N - 1
            DO 10 I = 0, M - 1
               B( I, J ) = ZERO
            } // 10
         } // 20
         RETURN
      }

      if ( LSIDE ) {

         // SIDE = 'L'

         // A is M-by-M.
         // If M is odd, set NISODD = .TRUE., and M1 and M2.
         // If M is even, NISODD = .FALSE., and M.

         if ( MOD( M, 2 ).EQ.0 ) {
            MISODD = .FALSE.
            K = M / 2
         } else {
            MISODD = .TRUE.
            if ( LOWER ) {
               M2 = M / 2
               M1 = M - M2
            } else {
               M1 = M / 2
               M2 = M - M1
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

                     if ( M.EQ.1 ) {
                        strsm('L', 'L', 'N', DIAG, M1, N, ALPHA, A, M, B, LDB );
                     } else {
                        strsm('L', 'L', 'N', DIAG, M1, N, ALPHA, A( 0 ), M, B, LDB )                         CALL SGEMM( 'N', 'N', M2, N, M1, -ONE, A( M1 ), M, B, LDB, ALPHA, B( M1, 0 ), LDB )                         CALL STRSM( 'L', 'U', 'T', DIAG, M2, N, ONE, A( M ), M, B( M1, 0 ), LDB );
                     }

                  } else {

                     // SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'L', and
                     // TRANS = 'T'

                     if ( M.EQ.1 ) {
                        strsm('L', 'L', 'T', DIAG, M1, N, ALPHA, A( 0 ), M, B, LDB );
                     } else {
                        strsm('L', 'U', 'N', DIAG, M2, N, ALPHA, A( M ), M, B( M1, 0 ), LDB )                         CALL SGEMM( 'T', 'N', M1, N, M2, -ONE, A( M1 ), M, B( M1, 0 ), LDB, ALPHA, B, LDB )                         CALL STRSM( 'L', 'L', 'T', DIAG, M1, N, ONE, A( 0 ), M, B, LDB );
                     }

                  }

               } else {

                  // SIDE  ='L', N is odd, TRANSR = 'N', and UPLO = 'U'

                  if ( .NOT.NOTRANS ) {

                     // SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'U', and
                     // TRANS = 'N'

                     strsm('L', 'L', 'N', DIAG, M1, N, ALPHA, A( M2 ), M, B, LDB )                      CALL SGEMM( 'T', 'N', M2, N, M1, -ONE, A( 0 ), M, B, LDB, ALPHA, B( M1, 0 ), LDB )                      CALL STRSM( 'L', 'U', 'T', DIAG, M2, N, ONE, A( M1 ), M, B( M1, 0 ), LDB );

                  } else {

                     // SIDE  ='L', N is odd, TRANSR = 'N', UPLO = 'U', and
                     // TRANS = 'T'

                     strsm('L', 'U', 'N', DIAG, M2, N, ALPHA, A( M1 ), M, B( M1, 0 ), LDB )                      CALL SGEMM( 'N', 'N', M1, N, M2, -ONE, A( 0 ), M, B( M1, 0 ), LDB, ALPHA, B, LDB )                      CALL STRSM( 'L', 'L', 'T', DIAG, M1, N, ONE, A( M2 ), M, B, LDB );

                  }

               }

            } else {

               // SIDE = 'L', N is odd, and TRANSR = 'T'

               if ( LOWER ) {

                  // SIDE  ='L', N is odd, TRANSR = 'T', and UPLO = 'L'

                  if ( NOTRANS ) {

                     // SIDE  ='L', N is odd, TRANSR = 'T', UPLO = 'L', and
                     // TRANS = 'N'

                     if ( M.EQ.1 ) {
                        strsm('L', 'U', 'T', DIAG, M1, N, ALPHA, A( 0 ), M1, B, LDB );
                     } else {
                        strsm('L', 'U', 'T', DIAG, M1, N, ALPHA, A( 0 ), M1, B, LDB )                         CALL SGEMM( 'T', 'N', M2, N, M1, -ONE, A( M1*M1 ), M1, B, LDB, ALPHA, B( M1, 0 ), LDB );
                        strsm('L', 'L', 'N', DIAG, M2, N, ONE, A( 1 ), M1, B( M1, 0 ), LDB );
                     }

                  } else {

                     // SIDE  ='L', N is odd, TRANSR = 'T', UPLO = 'L', and
                     // TRANS = 'T'

                     if ( M.EQ.1 ) {
                        strsm('L', 'U', 'N', DIAG, M1, N, ALPHA, A( 0 ), M1, B, LDB );
                     } else {
                        strsm('L', 'L', 'T', DIAG, M2, N, ALPHA, A( 1 ), M1, B( M1, 0 ), LDB )                         CALL SGEMM( 'N', 'N', M1, N, M2, -ONE, A( M1*M1 ), M1, B( M1, 0 ), LDB, ALPHA, B, LDB );
                        strsm('L', 'U', 'N', DIAG, M1, N, ONE, A( 0 ), M1, B, LDB );
                     }

                  }

               } else {

                  // SIDE  ='L', N is odd, TRANSR = 'T', and UPLO = 'U'

                  if ( .NOT.NOTRANS ) {

                     // SIDE  ='L', N is odd, TRANSR = 'T', UPLO = 'U', and
                     // TRANS = 'N'

                     strsm('L', 'U', 'T', DIAG, M1, N, ALPHA, A( M2*M2 ), M2, B, LDB )                      CALL SGEMM( 'N', 'N', M2, N, M1, -ONE, A( 0 ), M2, B, LDB, ALPHA, B( M1, 0 ), LDB )                      CALL STRSM( 'L', 'L', 'N', DIAG, M2, N, ONE, A( M1*M2 ), M2, B( M1, 0 ), LDB );

                  } else {

                     // SIDE  ='L', N is odd, TRANSR = 'T', UPLO = 'U', and
                     // TRANS = 'T'

                     strsm('L', 'L', 'T', DIAG, M2, N, ALPHA, A( M1*M2 ), M2, B( M1, 0 ), LDB )                      CALL SGEMM( 'T', 'N', M1, N, M2, -ONE, A( 0 ), M2, B( M1, 0 ), LDB, ALPHA, B, LDB )                      CALL STRSM( 'L', 'U', 'N', DIAG, M1, N, ONE, A( M2*M2 ), M2, B, LDB );

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

                     strsm('L', 'L', 'N', DIAG, K, N, ALPHA, A( 1 ), M+1, B, LDB )                      CALL SGEMM( 'N', 'N', K, N, K, -ONE, A( K+1 ), M+1, B, LDB, ALPHA, B( K, 0 ), LDB )                      CALL STRSM( 'L', 'U', 'T', DIAG, K, N, ONE, A( 0 ), M+1, B( K, 0 ), LDB );

                  } else {

                     // SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'L',
                     // and TRANS = 'T'

                     strsm('L', 'U', 'N', DIAG, K, N, ALPHA, A( 0 ), M+1, B( K, 0 ), LDB )                      CALL SGEMM( 'T', 'N', K, N, K, -ONE, A( K+1 ), M+1, B( K, 0 ), LDB, ALPHA, B, LDB )                      CALL STRSM( 'L', 'L', 'T', DIAG, K, N, ONE, A( 1 ), M+1, B, LDB );

                  }

               } else {

                  // SIDE  ='L', N is even, TRANSR = 'N', and UPLO = 'U'

                  if ( .NOT.NOTRANS ) {

                     // SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'U',
                     // and TRANS = 'N'

                     strsm('L', 'L', 'N', DIAG, K, N, ALPHA, A( K+1 ), M+1, B, LDB )                      CALL SGEMM( 'T', 'N', K, N, K, -ONE, A( 0 ), M+1, B, LDB, ALPHA, B( K, 0 ), LDB )                      CALL STRSM( 'L', 'U', 'T', DIAG, K, N, ONE, A( K ), M+1, B( K, 0 ), LDB );

                  } else {

                     // SIDE  ='L', N is even, TRANSR = 'N', UPLO = 'U',
                     // and TRANS = 'T'
                     strsm('L', 'U', 'N', DIAG, K, N, ALPHA, A( K ), M+1, B( K, 0 ), LDB )                      CALL SGEMM( 'N', 'N', K, N, K, -ONE, A( 0 ), M+1, B( K, 0 ), LDB, ALPHA, B, LDB )                      CALL STRSM( 'L', 'L', 'T', DIAG, K, N, ONE, A( K+1 ), M+1, B, LDB );

                  }

               }

            } else {

               // SIDE = 'L', N is even, and TRANSR = 'T'

               if ( LOWER ) {

                  // SIDE  ='L', N is even, TRANSR = 'T', and UPLO = 'L'

                  if ( NOTRANS ) {

                     // SIDE  ='L', N is even, TRANSR = 'T', UPLO = 'L',
                     // and TRANS = 'N'

                     strsm('L', 'U', 'T', DIAG, K, N, ALPHA, A( K ), K, B, LDB )                      CALL SGEMM( 'T', 'N', K, N, K, -ONE, A( K*( K+1 ) ), K, B, LDB, ALPHA, B( K, 0 ), LDB );
                     strsm('L', 'L', 'N', DIAG, K, N, ONE, A( 0 ), K, B( K, 0 ), LDB );

                  } else {

                     // SIDE  ='L', N is even, TRANSR = 'T', UPLO = 'L',
                     // and TRANS = 'T'

                     strsm('L', 'L', 'T', DIAG, K, N, ALPHA, A( 0 ), K, B( K, 0 ), LDB )                      CALL SGEMM( 'N', 'N', K, N, K, -ONE, A( K*( K+1 ) ), K, B( K, 0 ), LDB, ALPHA, B, LDB );
                     strsm('L', 'U', 'N', DIAG, K, N, ONE, A( K ), K, B, LDB );

                  }

               } else {

                  // SIDE  ='L', N is even, TRANSR = 'T', and UPLO = 'U'

                  if ( .NOT.NOTRANS ) {

                     // SIDE  ='L', N is even, TRANSR = 'T', UPLO = 'U',
                     // and TRANS = 'N'

                     strsm('L', 'U', 'T', DIAG, K, N, ALPHA, A( K*( K+1 ) ), K, B, LDB )                      CALL SGEMM( 'N', 'N', K, N, K, -ONE, A( 0 ), K, B, LDB, ALPHA, B( K, 0 ), LDB )                      CALL STRSM( 'L', 'L', 'N', DIAG, K, N, ONE, A( K*K ), K, B( K, 0 ), LDB );

                  } else {

                     // SIDE  ='L', N is even, TRANSR = 'T', UPLO = 'U',
                     // and TRANS = 'T'

                     strsm('L', 'L', 'T', DIAG, K, N, ALPHA, A( K*K ), K, B( K, 0 ), LDB )                      CALL SGEMM( 'T', 'N', K, N, K, -ONE, A( 0 ), K, B( K, 0 ), LDB, ALPHA, B, LDB )                      CALL STRSM( 'L', 'U', 'N', DIAG, K, N, ONE, A( K*( K+1 ) ), K, B, LDB );

                  }

               }

            }

         }

      } else {

         // SIDE = 'R'

         // A is N-by-N.
         // If N is odd, set NISODD = .TRUE., and N1 and N2.
         // If N is even, NISODD = .FALSE., and K.

         if ( MOD( N, 2 ).EQ.0 ) {
            NISODD = .FALSE.
            K = N / 2
         } else {
            NISODD = .TRUE.
            if ( LOWER ) {
               N2 = N / 2
               N1 = N - N2
            } else {
               N1 = N / 2
               N2 = N - N1
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

                     strsm('R', 'U', 'T', DIAG, M, N2, ALPHA, A( N ), N, B( 0, N1 ), LDB )                      CALL SGEMM( 'N', 'N', M, N1, N2, -ONE, B( 0, N1 ), LDB, A( N1 ), N, ALPHA, B( 0, 0 ), LDB );
                     strsm('R', 'L', 'N', DIAG, M, N1, ONE, A( 0 ), N, B( 0, 0 ), LDB );

                  } else {

                     // SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'L', and
                     // TRANS = 'T'

                     strsm('R', 'L', 'T', DIAG, M, N1, ALPHA, A( 0 ), N, B( 0, 0 ), LDB )                      CALL SGEMM( 'N', 'T', M, N2, N1, -ONE, B( 0, 0 ), LDB, A( N1 ), N, ALPHA, B( 0, N1 ), LDB );
                     strsm('R', 'U', 'N', DIAG, M, N2, ONE, A( N ), N, B( 0, N1 ), LDB );

                  }

               } else {

                  // SIDE  ='R', N is odd, TRANSR = 'N', and UPLO = 'U'

                  if ( NOTRANS ) {

                     // SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'U', and
                     // TRANS = 'N'

                     strsm('R', 'L', 'T', DIAG, M, N1, ALPHA, A( N2 ), N, B( 0, 0 ), LDB )                      CALL SGEMM( 'N', 'N', M, N2, N1, -ONE, B( 0, 0 ), LDB, A( 0 ), N, ALPHA, B( 0, N1 ), LDB );
                     strsm('R', 'U', 'N', DIAG, M, N2, ONE, A( N1 ), N, B( 0, N1 ), LDB );

                  } else {

                     // SIDE  ='R', N is odd, TRANSR = 'N', UPLO = 'U', and
                     // TRANS = 'T'

                     strsm('R', 'U', 'T', DIAG, M, N2, ALPHA, A( N1 ), N, B( 0, N1 ), LDB )                      CALL SGEMM( 'N', 'T', M, N1, N2, -ONE, B( 0, N1 ), LDB, A( 0 ), N, ALPHA, B( 0, 0 ), LDB )                      CALL STRSM( 'R', 'L', 'N', DIAG, M, N1, ONE, A( N2 ), N, B( 0, 0 ), LDB );

                  }

               }

            } else {

               // SIDE = 'R', N is odd, and TRANSR = 'T'

               if ( LOWER ) {

                  // SIDE  ='R', N is odd, TRANSR = 'T', and UPLO = 'L'

                  if ( NOTRANS ) {

                     // SIDE  ='R', N is odd, TRANSR = 'T', UPLO = 'L', and
                     // TRANS = 'N'

                     strsm('R', 'L', 'N', DIAG, M, N2, ALPHA, A( 1 ), N1, B( 0, N1 ), LDB )                      CALL SGEMM( 'N', 'T', M, N1, N2, -ONE, B( 0, N1 ), LDB, A( N1*N1 ), N1, ALPHA, B( 0, 0 ), LDB );
                     strsm('R', 'U', 'T', DIAG, M, N1, ONE, A( 0 ), N1, B( 0, 0 ), LDB );

                  } else {

                     // SIDE  ='R', N is odd, TRANSR = 'T', UPLO = 'L', and
                     // TRANS = 'T'

                     strsm('R', 'U', 'N', DIAG, M, N1, ALPHA, A( 0 ), N1, B( 0, 0 ), LDB )                      CALL SGEMM( 'N', 'N', M, N2, N1, -ONE, B( 0, 0 ), LDB, A( N1*N1 ), N1, ALPHA, B( 0, N1 ), LDB );
                     strsm('R', 'L', 'T', DIAG, M, N2, ONE, A( 1 ), N1, B( 0, N1 ), LDB );

                  }

               } else {

                  // SIDE  ='R', N is odd, TRANSR = 'T', and UPLO = 'U'

                  if ( NOTRANS ) {

                     // SIDE  ='R', N is odd, TRANSR = 'T', UPLO = 'U', and
                     // TRANS = 'N'

                     strsm('R', 'U', 'N', DIAG, M, N1, ALPHA, A( N2*N2 ), N2, B( 0, 0 ), LDB )                      CALL SGEMM( 'N', 'T', M, N2, N1, -ONE, B( 0, 0 ), LDB, A( 0 ), N2, ALPHA, B( 0, N1 ), LDB );
                     strsm('R', 'L', 'T', DIAG, M, N2, ONE, A( N1*N2 ), N2, B( 0, N1 ), LDB );

                  } else {

                     // SIDE  ='R', N is odd, TRANSR = 'T', UPLO = 'U', and
                     // TRANS = 'T'

                     strsm('R', 'L', 'N', DIAG, M, N2, ALPHA, A( N1*N2 ), N2, B( 0, N1 ), LDB )                      CALL SGEMM( 'N', 'N', M, N1, N2, -ONE, B( 0, N1 ), LDB, A( 0 ), N2, ALPHA, B( 0, 0 ), LDB );
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

                     strsm('R', 'U', 'T', DIAG, M, K, ALPHA, A( 0 ), N+1, B( 0, K ), LDB )                      CALL SGEMM( 'N', 'N', M, K, K, -ONE, B( 0, K ), LDB, A( K+1 ), N+1, ALPHA, B( 0, 0 ), LDB );
                     strsm('R', 'L', 'N', DIAG, M, K, ONE, A( 1 ), N+1, B( 0, 0 ), LDB );

                  } else {

                     // SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'L',
                     // and TRANS = 'T'

                     strsm('R', 'L', 'T', DIAG, M, K, ALPHA, A( 1 ), N+1, B( 0, 0 ), LDB )                      CALL SGEMM( 'N', 'T', M, K, K, -ONE, B( 0, 0 ), LDB, A( K+1 ), N+1, ALPHA, B( 0, K ), LDB );
                     strsm('R', 'U', 'N', DIAG, M, K, ONE, A( 0 ), N+1, B( 0, K ), LDB );

                  }

               } else {

                  // SIDE  ='R', N is even, TRANSR = 'N', and UPLO = 'U'

                  if ( NOTRANS ) {

                     // SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'U',
                     // and TRANS = 'N'

                     strsm('R', 'L', 'T', DIAG, M, K, ALPHA, A( K+1 ), N+1, B( 0, 0 ), LDB )                      CALL SGEMM( 'N', 'N', M, K, K, -ONE, B( 0, 0 ), LDB, A( 0 ), N+1, ALPHA, B( 0, K ), LDB );
                     strsm('R', 'U', 'N', DIAG, M, K, ONE, A( K ), N+1, B( 0, K ), LDB );

                  } else {

                     // SIDE  ='R', N is even, TRANSR = 'N', UPLO = 'U',
                     // and TRANS = 'T'

                     strsm('R', 'U', 'T', DIAG, M, K, ALPHA, A( K ), N+1, B( 0, K ), LDB )                      CALL SGEMM( 'N', 'T', M, K, K, -ONE, B( 0, K ), LDB, A( 0 ), N+1, ALPHA, B( 0, 0 ), LDB );
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

                     strsm('R', 'L', 'N', DIAG, M, K, ALPHA, A( 0 ), K, B( 0, K ), LDB )                      CALL SGEMM( 'N', 'T', M, K, K, -ONE, B( 0, K ), LDB, A( ( K+1 )*K ), K, ALPHA, B( 0, 0 ), LDB );
                     strsm('R', 'U', 'T', DIAG, M, K, ONE, A( K ), K, B( 0, 0 ), LDB );

                  } else {

                     // SIDE  ='R', N is even, TRANSR = 'T', UPLO = 'L',
                     // and TRANS = 'T'

                     strsm('R', 'U', 'N', DIAG, M, K, ALPHA, A( K ), K, B( 0, 0 ), LDB )                      CALL SGEMM( 'N', 'N', M, K, K, -ONE, B( 0, 0 ), LDB, A( ( K+1 )*K ), K, ALPHA, B( 0, K ), LDB );
                     strsm('R', 'L', 'T', DIAG, M, K, ONE, A( 0 ), K, B( 0, K ), LDB );

                  }

               } else {

                  // SIDE  ='R', N is even, TRANSR = 'T', and UPLO = 'U'

                  if ( NOTRANS ) {

                     // SIDE  ='R', N is even, TRANSR = 'T', UPLO = 'U',
                     // and TRANS = 'N'

                     strsm('R', 'U', 'N', DIAG, M, K, ALPHA, A( ( K+1 )*K ), K, B( 0, 0 ), LDB )                      CALL SGEMM( 'N', 'T', M, K, K, -ONE, B( 0, 0 ), LDB, A( 0 ), K, ALPHA, B( 0, K ), LDB )                      CALL STRSM( 'R', 'L', 'T', DIAG, M, K, ONE, A( K*K ), K, B( 0, K ), LDB );

                  } else {

                     // SIDE  ='R', N is even, TRANSR = 'T', UPLO = 'U',
                     // and TRANS = 'T'

                     strsm('R', 'L', 'N', DIAG, M, K, ALPHA, A( K*K ), K, B( 0, K ), LDB )                      CALL SGEMM( 'N', 'N', M, K, K, -ONE, B( 0, K ), LDB, A( 0 ), K, ALPHA, B( 0, 0 ), LDB )                      CALL STRSM( 'R', 'U', 'T', DIAG, M, K, ONE, A( ( K+1 )*K ), K, B( 0, 0 ), LDB );

                  }

               }

            }

         }
      }

      RETURN

      // End of STFSM

      }
