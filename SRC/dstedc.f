      SUBROUTINE DSTEDC( COMPZ, N, D, E, Z, LDZ, WORK, LWORK, IWORK, LIWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             COMPZ;
      int                INFO, LDZ, LIWORK, LWORK, N;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             D( * ), E( * ), WORK( * ), Z( LDZ, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TWO;
      const              ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                FINISH, I, ICOMPZ, II, J, K, LGN, LIWMIN, LWMIN, M, SMLSIZ, START, STOREZ, STRTRW;
      double             EPS, ORGNRM, P, TINY;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      double             DLAMCH, DLANST;
      // EXTERNAL LSAME, ILAENV, DLAMCH, DLANST
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGEMM, DLACPY, DLAED0, DLASCL, DLASET, DLASRT, DSTEQR, DSTERF, DSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, INT, LOG, MAX, MOD, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )

      if ( LSAME( COMPZ, 'N' ) ) {
         ICOMPZ = 0
      } else if ( LSAME( COMPZ, 'V' ) ) {
         ICOMPZ = 1
      } else if ( LSAME( COMPZ, 'I' ) ) {
         ICOMPZ = 2
      } else {
         ICOMPZ = -1
      }
      if ( ICOMPZ.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1, N ) ) ) {
         INFO = -6
      }

      if ( INFO.EQ.0 ) {

         // Compute the workspace requirements

         SMLSIZ = ILAENV( 9, 'DSTEDC', ' ', 0, 0, 0, 0 )
         if ( N.LE.1 .OR. ICOMPZ.EQ.0 ) {
            LIWMIN = 1
            LWMIN = 1
         } else if ( N.LE.SMLSIZ ) {
            LIWMIN = 1
            LWMIN = 2*( N - 1 )
         } else {
            LGN = INT( LOG( DBLE( N ) )/LOG( TWO ) )
            if (2**LGN.LT.N) LGN = LGN + 1             IF( 2**LGN.LT.N ) LGN = LGN + 1;
            if ( ICOMPZ.EQ.1 ) {
               LWMIN = 1 + 3*N + 2*N*LGN + 4*N**2
               LIWMIN = 6 + 6*N + 5*N*LGN
            } else if ( ICOMPZ.EQ.2 ) {
               LWMIN = 1 + 4*N + N**2
               LIWMIN = 3 + 5*N
            }
         }
         WORK( 1 ) = LWMIN
         IWORK( 1 ) = LIWMIN

         if ( LWORK.LT.LWMIN .AND. .NOT. LQUERY ) {
            INFO = -8
         } else if ( LIWORK.LT.LIWMIN .AND. .NOT. LQUERY ) {
            INFO = -10
         }
      }

      if ( INFO.NE.0 ) {
         xerbla('DSTEDC', -INFO );
         RETURN
      } else if (LQUERY) {
         RETURN
      }

      // Quick return if possible

      if (N.EQ.0) RETURN;
      if ( N.EQ.1 ) {
         if (ICOMPZ.NE.0) Z( 1, 1 ) = ONE;
         RETURN
      }

      // If the following conditional clause is removed, then the routine
      // will use the Divide and Conquer routine to compute only the
      // eigenvalues, which requires (3N + 3N**2) real workspace and
      // (2 + 5N + 2N lg(N)) integer workspace.
      // Since on many architectures DSTERF is much faster than any other
      // algorithm for finding eigenvalues only, it is used here
      // as the default. If the conditional clause is removed, then
      // information on the size of workspace needs to be changed.

      // If COMPZ = 'N', use DSTERF to compute the eigenvalues.

      if ( ICOMPZ.EQ.0 ) {
         dsterf(N, D, E, INFO );
         GO TO 50
      }

      // If N is smaller than the minimum divide size (SMLSIZ+1), then
      // solve the problem with another solver.

      if ( N.LE.SMLSIZ ) {

         dsteqr(COMPZ, N, D, E, Z, LDZ, WORK, INFO );

      } else {

         // If COMPZ = 'V', the Z matrix must be stored elsewhere for later
         // use.

         if ( ICOMPZ.EQ.1 ) {
            STOREZ = 1 + N*N
         } else {
            STOREZ = 1
         }

         if ( ICOMPZ.EQ.2 ) {
            dlaset('Full', N, N, ZERO, ONE, Z, LDZ );
         }

         // Scale.

         ORGNRM = DLANST( 'M', N, D, E )
         if (ORGNRM.EQ.ZERO) GO TO 50;

         EPS = DLAMCH( 'Epsilon' )

         START = 1

         // while ( START <= N )

         } // 10
         if ( START.LE.N ) {

            // Let FINISH be the position of the next subdiagonal entry
            // such that E( FINISH ) <= TINY or FINISH = N if no such
            // subdiagonal exists.  The matrix identified by the elements
            // between START and FINISH constitutes an independent
            // sub-problem.

            FINISH = START
            } // 20
            if ( FINISH.LT.N ) {
               TINY = EPS*SQRT( ABS( D( FINISH ) ) )* SQRT( ABS( D( FINISH+1 ) ) )
               if ( ABS( E( FINISH ) ).GT.TINY ) {
                  FINISH = FINISH + 1
                  GO TO 20
               }
            }

            // (Sub) Problem determined.  Compute its size and solve it.

            M = FINISH - START + 1
            if ( M.EQ.1 ) {
               START = FINISH + 1
               GO TO 10
            }
            if ( M.GT.SMLSIZ ) {

               // Scale.

               ORGNRM = DLANST( 'M', M, D( START ), E( START ) )
               dlascl('G', 0, 0, ORGNRM, ONE, M, 1, D( START ), M, INFO );
               dlascl('G', 0, 0, ORGNRM, ONE, M-1, 1, E( START ), M-1, INFO );

               if ( ICOMPZ.EQ.1 ) {
                  STRTRW = 1
               } else {
                  STRTRW = START
               }
               dlaed0(ICOMPZ, N, M, D( START ), E( START ), Z( STRTRW, START ), LDZ, WORK( 1 ), N, WORK( STOREZ ), IWORK, INFO );
               if ( INFO.NE.0 ) {
                  INFO = ( INFO / ( M+1 )+START-1 )*( N+1 ) + MOD( INFO, ( M+1 ) ) + START - 1
                  GO TO 50
               }

               // Scale back.

               dlascl('G', 0, 0, ONE, ORGNRM, M, 1, D( START ), M, INFO );

            } else {
               if ( ICOMPZ.EQ.1 ) {

                  // Since QR won't update a Z matrix which is larger than
                  // the length of D, we must solve the sub-problem in a
                  // workspace and then multiply back into Z.

                  dsteqr('I', M, D( START ), E( START ), WORK, M, WORK( M*M+1 ), INFO );
                  dlacpy('A', N, M, Z( 1, START ), LDZ, WORK( STOREZ ), N );
                  dgemm('N', 'N', N, M, M, ONE, WORK( STOREZ ), N, WORK, M, ZERO, Z( 1, START ), LDZ );
               } else if ( ICOMPZ.EQ.2 ) {
                  dsteqr('I', M, D( START ), E( START ), Z( START, START ), LDZ, WORK, INFO );
               } else {
                  dsterf(M, D( START ), E( START ), INFO );
               }
               if ( INFO.NE.0 ) {
                  INFO = START*( N+1 ) + FINISH
                  GO TO 50
               }
            }

            START = FINISH + 1
            GO TO 10
         }

         // endwhile

         if ( ICOMPZ.EQ.0 ) {

           // Use Quick Sort

           dlasrt('I', N, D, INFO );

         } else {

           // Use Selection Sort to minimize swaps of eigenvectors

           for (II = 2; II <= N; II++) { // 40
              I = II - 1
              K = I
              P = D( I )
              for (J = II; J <= N; J++) { // 30
                 if ( D( J ).LT.P ) {
                    K = J
                    P = D( J )
                 }
              } // 30
              if ( K.NE.I ) {
                 D( K ) = D( I )
                 D( I ) = P
                 dswap(N, Z( 1, I ), 1, Z( 1, K ), 1 );
              }
           } // 40
         }
      }

      } // 50
      WORK( 1 ) = LWMIN
      IWORK( 1 ) = LIWMIN

      RETURN

      // End of DSTEDC

      }
