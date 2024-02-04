      void zstedc(COMPZ, N, D, E, Z, LDZ, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             COMPZ;
      int                INFO, LDZ, LIWORK, LRWORK, LWORK, N;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double             D( * ), E( * ), RWORK( * );
      Complex         WORK( * ), Z( LDZ, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TWO;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                FINISH, I, ICOMPZ, II, J, K, LGN, LIWMIN, LL, LRWMIN, LWMIN, M, SMLSIZ, START;
      double             EPS, ORGNRM, P, TINY;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV;
      //- double             DLAMCH, DLANST;
      // EXTERNAL lsame, ILAENV, DLAMCH, DLANST
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASCL, DLASET, DSTEDC, DSTEQR, DSTERF, XERBLA, ZLACPY, ZLACRM, ZLAED0, ZSTEQR, ZSWAP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, INT, LOG, MAX, MOD, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      LQUERY = ( LWORK == -1 || LRWORK == -1 || LIWORK == -1 );

      if ( lsame( COMPZ, 'N' ) ) {
         ICOMPZ = 0;
      } else if ( lsame( COMPZ, 'V' ) ) {
         ICOMPZ = 1;
      } else if ( lsame( COMPZ, 'I' ) ) {
         ICOMPZ = 2;
      } else {
         ICOMPZ = -1;
      }
      if ( ICOMPZ < 0 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( ( LDZ < 1 ) || ( ICOMPZ > 0 && LDZ < max( 1, N ) ) ) {
         INFO = -6;
      }

      if ( INFO == 0 ) {

         // Compute the workspace requirements

         SMLSIZ = ILAENV( 9, 'ZSTEDC', ' ', 0, 0, 0, 0 );
         if ( N <= 1 || ICOMPZ == 0 ) {
            LWMIN = 1;
            LIWMIN = 1;
            LRWMIN = 1;
         } else if ( N <= SMLSIZ ) {
            LWMIN = 1;
            LIWMIN = 1;
            LRWMIN = 2*( N - 1 );
         } else if ( ICOMPZ == 1 ) {
            LGN = INT( LOG( N.toDouble() ) / LOG( TWO ) );
            if (2**LGN < N) LGN = LGN + 1;
            IF( 2**LGN < N ) LGN = LGN + 1;
            LWMIN = N*N;
            LRWMIN = 1 + 3*N + 2*N*LGN + 4*N**2;
            LIWMIN = 6 + 6*N + 5*N*LGN;
         } else if ( ICOMPZ == 2 ) {
            LWMIN = 1;
            LRWMIN = 1 + 4*N + 2*N**2;
            LIWMIN = 3 + 5*N;
         }
         WORK[1] = LWMIN;
         RWORK[1] = LRWMIN;
         IWORK[1] = LIWMIN;

         if ( LWORK < LWMIN && !LQUERY ) {
            INFO = -8;
         } else if ( LRWORK < LRWMIN && !LQUERY ) {
            INFO = -10;
         } else if ( LIWORK < LIWMIN && !LQUERY ) {
            INFO = -12;
         }
      }

      if ( INFO != 0 ) {
         xerbla('ZSTEDC', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;
      if ( N == 1 ) {
         if (ICOMPZ != 0) Z( 1, 1 ) = ONE;
         return;
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

      if ( ICOMPZ == 0 ) {
         dsterf(N, D, E, INFO );
         GO TO 70;
      }

      // If N is smaller than the minimum divide size (SMLSIZ+1), then
      // solve the problem with another solver.

      if ( N <= SMLSIZ ) {

         zsteqr(COMPZ, N, D, E, Z, LDZ, RWORK, INFO );

      } else {

         // If COMPZ = 'I', we simply call DSTEDC instead.

         if ( ICOMPZ == 2 ) {
            dlaset('Full', N, N, ZERO, ONE, RWORK, N );
            LL = N*N + 1;
            dstedc('I', N, D, E, RWORK, N, RWORK( LL ), LRWORK-LL+1, IWORK, LIWORK, INFO );
            for (J = 1; J <= N; J++) { // 20
               for (I = 1; I <= N; I++) { // 10
                  Z[I, J] = RWORK( ( J-1 )*N+I );
               } // 10
            } // 20
            GO TO 70;
         }

         // From now on, only option left to be handled is COMPZ = 'V',
         // i.e. ICOMPZ = 1.

         // Scale.

         ORGNRM = DLANST( 'M', N, D, E );
         if (ORGNRM == ZERO) GO TO 70;

         EPS = DLAMCH( 'Epsilon' );

         START = 1;

         // while ( START <= N )

         } // 30
         if ( START <= N ) {

            // Let FINISH be the position of the next subdiagonal entry
            // such that E( FINISH ) <= TINY or FINISH = N if no such
            // subdiagonal exists.  The matrix identified by the elements
            // between START and FINISH constitutes an independent
            // sub-problem.

            FINISH = START;
            } // 40
            if ( FINISH < N ) {
               TINY = EPS*sqrt( ( D( FINISH ) ) ).abs()* sqrt( ( D( FINISH+1 ) ) ).abs();
               if ( ( E( FINISH ) ).abs() > TINY ) {
                  FINISH = FINISH + 1;
                  GO TO 40;
               }
            }

            // (Sub) Problem determined.  Compute its size and solve it.

            M = FINISH - START + 1;
            if ( M > SMLSIZ ) {

               // Scale.

               ORGNRM = DLANST( 'M', M, D( START ), E( START ) );
               dlascl('G', 0, 0, ORGNRM, ONE, M, 1, D( START ), M, INFO );
               dlascl('G', 0, 0, ORGNRM, ONE, M-1, 1, E( START ), M-1, INFO );

               zlaed0(N, M, D( START ), E( START ), Z( 1, START ), LDZ, WORK, N, RWORK, IWORK, INFO );
               if ( INFO > 0 ) {
                  INFO = ( INFO / ( M+1 )+START-1 )*( N+1 ) + (INFO % ( M+1) ) + START - 1;
                  GO TO 70;
               }

               // Scale back.

               dlascl('G', 0, 0, ONE, ORGNRM, M, 1, D( START ), M, INFO );

            } else {
               dsteqr('I', M, D( START ), E( START ), RWORK, M, RWORK( M*M+1 ), INFO );
               zlacrm(N, M, Z( 1, START ), LDZ, RWORK, M, WORK, N, RWORK( M*M+1 ) );
               zlacpy('A', N, M, WORK, N, Z( 1, START ), LDZ );
               if ( INFO > 0 ) {
                  INFO = START*( N+1 ) + FINISH;
                  GO TO 70;
               }
            }

            START = FINISH + 1;
            GO TO 30;
         }

         // endwhile


         // Use Selection Sort to minimize swaps of eigenvectors

         for (II = 2; II <= N; II++) { // 60
           I = II - 1;
           K = I;
           P = D( I );
           for (J = II; J <= N; J++) { // 50
              if ( D( J ) < P ) {
                 K = J;
                 P = D( J );
              }
           } // 50
           if ( K != I ) {
              D[K] = D( I );
              D[I] = P;
              zswap(N, Z( 1, I ), 1, Z( 1, K ), 1 );
           }
         } // 60
      }

      } // 70
      WORK[1] = LWMIN;
      RWORK[1] = LRWMIN;
      IWORK[1] = LIWMIN;

      return;
      }