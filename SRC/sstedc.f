      SUBROUTINE SSTEDC( COMPZ, N, D, E, Z, LDZ, WORK, LWORK, IWORK, LIWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             COMPZ;
      int                INFO, LDZ, LIWORK, LWORK, N;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL               D( * ), E( * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TWO
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                FINISH, I, ICOMPZ, II, J, K, LGN, LIWMIN, LWMIN, M, SMLSIZ, START, STOREZ, STRTRW;
      REAL               EPS, ORGNRM, P, TINY
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               SLAMCH, SLANST, SROUNDUP_LWORK
      // EXTERNAL ILAENV, LSAME, SLAMCH, SLANST, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SLACPY, SLAED0, SLASCL, SLASET, SLASRT, SSTEQR, SSTERF, SSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, INT, LOG, MAX, MOD, REAL, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      LQUERY = ( LWORK == -1 || LIWORK == -1 )

      if ( LSAME( COMPZ, 'N' ) ) {
         ICOMPZ = 0
      } else if ( LSAME( COMPZ, 'V' ) ) {
         ICOMPZ = 1
      } else if ( LSAME( COMPZ, 'I' ) ) {
         ICOMPZ = 2
      } else {
         ICOMPZ = -1
      }
      if ( ICOMPZ < 0 ) {
         INFO = -1
      } else if ( N < 0 ) {
         INFO = -2
      } else if ( ( LDZ < 1 ) || ( ICOMPZ > 0 && LDZ < MAX( 1, N ) ) ) {
         INFO = -6
      }

      if ( INFO == 0 ) {

         // Compute the workspace requirements

         SMLSIZ = ILAENV( 9, 'SSTEDC', ' ', 0, 0, 0, 0 )
         if ( N <= 1 || ICOMPZ == 0 ) {
            LIWMIN = 1
            LWMIN = 1
         } else if ( N <= SMLSIZ ) {
            LIWMIN = 1
            LWMIN = 2*( N - 1 )
         } else {
            LGN = INT( LOG( REAL( N ) )/LOG( TWO ) )
            if (2**LGN < N) LGN = LGN + 1             IF( 2**LGN < N ) LGN = LGN + 1;
            if ( ICOMPZ == 1 ) {
               LWMIN = 1 + 3*N + 2*N*LGN + 4*N**2
               LIWMIN = 6 + 6*N + 5*N*LGN
            } else if ( ICOMPZ == 2 ) {
               LWMIN = 1 + 4*N + N**2
               LIWMIN = 3 + 5*N
            }
         }
         WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
         IWORK( 1 ) = LIWMIN

         if ( LWORK < LWMIN && !LQUERY ) {
            INFO = -8
         } else if ( LIWORK < LIWMIN && !LQUERY ) {
            INFO = -10
         }
      }

      if ( INFO != 0 ) {
         xerbla('SSTEDC', -INFO );
         RETURN
      } else if (LQUERY) {
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;
      if ( N == 1 ) {
         if (ICOMPZ != 0) Z( 1, 1 ) = ONE;
         RETURN
      }

      // If the following conditional clause is removed, then the routine
      // will use the Divide and Conquer routine to compute only the
      // eigenvalues, which requires (3N + 3N**2) real workspace and
      // (2 + 5N + 2N lg(N)) integer workspace.
      // Since on many architectures SSTERF is much faster than any other
      // algorithm for finding eigenvalues only, it is used here
      // as the default. If the conditional clause is removed, then
      // information on the size of workspace needs to be changed.

      // If COMPZ = 'N', use SSTERF to compute the eigenvalues.

      if ( ICOMPZ == 0 ) {
         ssterf(N, D, E, INFO );
         GO TO 50
      }

      // If N is smaller than the minimum divide size (SMLSIZ+1), then
      // solve the problem with another solver.

      if ( N <= SMLSIZ ) {

         ssteqr(COMPZ, N, D, E, Z, LDZ, WORK, INFO );

      } else {

         // If COMPZ = 'V', the Z matrix must be stored elsewhere for later
         // use.

         if ( ICOMPZ == 1 ) {
            STOREZ = 1 + N*N
         } else {
            STOREZ = 1
         }

         if ( ICOMPZ == 2 ) {
            slaset('Full', N, N, ZERO, ONE, Z, LDZ );
         }

         // Scale.

         ORGNRM = SLANST( 'M', N, D, E )
         if (ORGNRM == ZERO) GO TO 50;

         EPS = SLAMCH( 'Epsilon' )

         START = 1

         // while ( START <= N )

         } // 10
         if ( START <= N ) {

            // Let FINISH be the position of the next subdiagonal entry
            // such that E( FINISH ) <= TINY or FINISH = N if no such
            // subdiagonal exists.  The matrix identified by the elements
            // between START and FINISH constitutes an independent
            // sub-problem.

            FINISH = START
            } // 20
            if ( FINISH < N ) {
               TINY = EPS*SQRT( ABS( D( FINISH ) ) )* SQRT( ABS( D( FINISH+1 ) ) )
               if ( ABS( E( FINISH ) ) > TINY ) {
                  FINISH = FINISH + 1
                  GO TO 20
               }
            }

            // (Sub) Problem determined.  Compute its size and solve it.

            M = FINISH - START + 1
            if ( M == 1 ) {
               START = FINISH + 1
               GO TO 10
            }
            if ( M > SMLSIZ ) {

               // Scale.

               ORGNRM = SLANST( 'M', M, D( START ), E( START ) )
               slascl('G', 0, 0, ORGNRM, ONE, M, 1, D( START ), M, INFO );
               slascl('G', 0, 0, ORGNRM, ONE, M-1, 1, E( START ), M-1, INFO );

               if ( ICOMPZ == 1 ) {
                  STRTRW = 1
               } else {
                  STRTRW = START
               }
               slaed0(ICOMPZ, N, M, D( START ), E( START ), Z( STRTRW, START ), LDZ, WORK( 1 ), N, WORK( STOREZ ), IWORK, INFO );
               if ( INFO != 0 ) {
                  INFO = ( INFO / ( M+1 )+START-1 )*( N+1 ) + MOD( INFO, ( M+1 ) ) + START - 1
                  GO TO 50
               }

               // Scale back.

               slascl('G', 0, 0, ONE, ORGNRM, M, 1, D( START ), M, INFO );

            } else {
               if ( ICOMPZ == 1 ) {

                  // Since QR won't update a Z matrix which is larger than
                  // the length of D, we must solve the sub-problem in a
                  // workspace and then multiply back into Z.

                  ssteqr('I', M, D( START ), E( START ), WORK, M, WORK( M*M+1 ), INFO );
                  slacpy('A', N, M, Z( 1, START ), LDZ, WORK( STOREZ ), N );
                  sgemm('N', 'N', N, M, M, ONE, WORK( STOREZ ), N, WORK, M, ZERO, Z( 1, START ), LDZ );
               } else if ( ICOMPZ == 2 ) {
                  ssteqr('I', M, D( START ), E( START ), Z( START, START ), LDZ, WORK, INFO );
               } else {
                  ssterf(M, D( START ), E( START ), INFO );
               }
               if ( INFO != 0 ) {
                  INFO = START*( N+1 ) + FINISH
                  GO TO 50
               }
            }

            START = FINISH + 1
            GO TO 10
         }

         // endwhile

         if ( ICOMPZ == 0 ) {

           // Use Quick Sort

           slasrt('I', N, D, INFO );

         } else {

           // Use Selection Sort to minimize swaps of eigenvectors

           for (II = 2; II <= N; II++) { // 40
              I = II - 1
              K = I
              P = D( I )
              for (J = II; J <= N; J++) { // 30
                 if ( D( J ) < P ) {
                    K = J
                    P = D( J )
                 }
              } // 30
              if ( K != I ) {
                 D( K ) = D( I )
                 D( I ) = P
                 sswap(N, Z( 1, I ), 1, Z( 1, K ), 1 );
              }
           } // 40
         }
      }

      } // 50
      WORK( 1 ) = SROUNDUP_LWORK(LWMIN)
      IWORK( 1 ) = LIWMIN

      RETURN

      // End of SSTEDC

      }
