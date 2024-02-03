      RECURSIVE SUBROUTINE CLAQZ0( WANTS, WANTQ, WANTZ, N, ILO, IHI, A, LDA, B, LDB, ALPHA, BETA, Q, LDQ, Z, LDZ, WORK, LWORK, RWORK, REC, INFO );
      // IMPLICIT NONE

      // Arguments
      String   , INTENT( IN ) :: WANTS, WANTQ, WANTZ;
      int    , INTENT( IN ) :: N, ILO, IHI, LDA, LDB, LDQ, LDZ, LWORK, REC;
      int    , INTENT( OUT ) :: INFO;
      COMPLEX, INTENT( INOUT ) :: A( LDA, * ), B( LDB, * ), Q( LDQ, * ), Z( LDZ, * ), ALPHA( * ), BETA( * ), WORK( * );
      REAL, INTENT( OUT ) :: RWORK( * );

      // Parameters
      COMPLEX         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      REAL :: ZERO, ONE, HALF;
      const    ZERO = 0.0, ONE = 1.0, HALF = 0.5 ;

      // Local scalars
      REAL :: SMLNUM, ULP, SAFMIN, SAFMAX, C1, TEMPR, BNORM, BTOL;
      COMPLEX :: ESHIFT, S1, TEMP;
      int     :: ISTART, ISTOP, IITER, MAXIT, ISTART2, K, LD, NSHIFTS, NBLOCK, NW, NMIN, NIBBLE, N_UNDEFLATED, N_DEFLATED, NS, SWEEP_INFO, SHIFTPOS, LWORKREQ, K2, ISTARTM, ISTOPM, IWANTS, IWANTQ, IWANTZ, NORM_INFO, AED_INFO, NWR, NBR, NSR, ITEMP1, ITEMP2, RCOST;
      bool    :: ILSCHUR, ILQ, ILZ;
      String    :: JBCMPZ*3;

      // External Functions
      // EXTERNAL :: XERBLA, CHGEQZ, CLAQZ2, CLAQZ3, CLASET, CLARTG, CROT
      REAL, EXTERNAL :: SLAMCH, CLANHS;
      bool   , EXTERNAL :: LSAME;
      int    , EXTERNAL :: ILAENV;


      // Decode wantS,wantQ,wantZ

      if ( LSAME( WANTS, 'E' ) ) {
         ILSCHUR = false;
         IWANTS = 1;
      } else if ( LSAME( WANTS, 'S' ) ) {
         ILSCHUR = true;
         IWANTS = 2;
      } else {
         IWANTS = 0;
      }

      if ( LSAME( WANTQ, 'N' ) ) {
         ILQ = false;
         IWANTQ = 1;
      } else if ( LSAME( WANTQ, 'V' ) ) {
         ILQ = true;
         IWANTQ = 2;
      } else if ( LSAME( WANTQ, 'I' ) ) {
         ILQ = true;
         IWANTQ = 3;
      } else {
         IWANTQ = 0;
      }

      if ( LSAME( WANTZ, 'N' ) ) {
         ILZ = false;
         IWANTZ = 1;
      } else if ( LSAME( WANTZ, 'V' ) ) {
         ILZ = true;
         IWANTZ = 2;
      } else if ( LSAME( WANTZ, 'I' ) ) {
         ILZ = true;
         IWANTZ = 3;
      } else {
         IWANTZ = 0;
      }

      // Check Argument Values

      INFO = 0;
      if ( IWANTS == 0 ) {
         INFO = -1;
      } else if ( IWANTQ == 0 ) {
         INFO = -2;
      } else if ( IWANTZ == 0 ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( ILO < 1 ) {
         INFO = -5;
      } else if ( IHI > N || IHI < ILO-1 ) {
         INFO = -6;
      } else if ( LDA < N ) {
         INFO = -8;
      } else if ( LDB < N ) {
         INFO = -10;
      } else if ( LDQ < 1 || ( ILQ && LDQ < N ) ) {
         INFO = -15;
      } else if ( LDZ < 1 || ( ILZ && LDZ < N ) ) {
         INFO = -17;
      }
      if ( INFO != 0 ) {
         xerbla('CLAQZ0', -INFO );
         return;
      }


      // Quick return if possible

      if ( N <= 0 ) {
         WORK( 1 ) = REAL( 1 );
         return;
      }


      // Get the parameters

      JBCMPZ( 1:1 ) = WANTS;
      JBCMPZ( 2:2 ) = WANTQ;
      JBCMPZ( 3:3 ) = WANTZ;

      NMIN = ILAENV( 12, 'CLAQZ0', JBCMPZ, N, ILO, IHI, LWORK );

      NWR = ILAENV( 13, 'CLAQZ0', JBCMPZ, N, ILO, IHI, LWORK );
      NWR = max( 2, NWR );
      NWR = min( IHI-ILO+1, ( N-1 ) / 3, NWR );

      NIBBLE = ILAENV( 14, 'CLAQZ0', JBCMPZ, N, ILO, IHI, LWORK );

      NSR = ILAENV( 15, 'CLAQZ0', JBCMPZ, N, ILO, IHI, LWORK );
      NSR = min( NSR, ( N+6 ) / 9, IHI-ILO );
      NSR = max( 2, NSR-MOD( NSR, 2 ) );

      RCOST = ILAENV( 17, 'CLAQZ0', JBCMPZ, N, ILO, IHI, LWORK );
      ITEMP1 = INT( NSR/sqrt( 1+2*NSR/( REAL( RCOST )/100*N ) ) );
      ITEMP1 = ( ( ITEMP1-1 )/4 )*4+4;
      NBR = NSR+ITEMP1;

      if ( N < NMIN || REC >= 2 ) {
         chgeqz(WANTS, WANTQ, WANTZ, N, ILO, IHI, A, LDA, B, LDB, ALPHA, BETA, Q, LDQ, Z, LDZ, WORK, LWORK, RWORK, INFO );
         return;
      }


      // Find out required workspace


      // Workspace query to CLAQZ2
      NW = max( NWR, NMIN );
      claqz2(ILSCHUR, ILQ, ILZ, N, ILO, IHI, NW, A, LDA, B, LDB, Q, LDQ, Z, LDZ, N_UNDEFLATED, N_DEFLATED, ALPHA, BETA, WORK, NW, WORK, NW, WORK, -1, RWORK, REC, AED_INFO );
      ITEMP1 = INT( WORK( 1 ) );
      // Workspace query to CLAQZ3
      claqz3(ILSCHUR, ILQ, ILZ, N, ILO, IHI, NSR, NBR, ALPHA, BETA, A, LDA, B, LDB, Q, LDQ, Z, LDZ, WORK, NBR, WORK, NBR, WORK, -1, SWEEP_INFO );
      ITEMP2 = INT( WORK( 1 ) );

      LWORKREQ = max( ITEMP1+2*NW**2, ITEMP2+2*NBR**2 );
      if ( LWORK == -1 ) {
         WORK( 1 ) = REAL( LWORKREQ );
         return;
      } else if ( LWORK < LWORKREQ ) {
         INFO = -19;
      }
      if ( INFO != 0 ) {
         xerbla('CLAQZ0', INFO );
         return;
      }

      // Initialize Q and Z

      if (IWANTQ == 3) claset( 'FULL', N, N, CZERO, CONE, Q, LDQ );
      IF( IWANTZ == 3 ) claset( 'FULL', N, N, CZERO, CONE, Z, LDZ );

      // Get machine constants
      SAFMIN = SLAMCH( 'SAFE MINIMUM' );
      SAFMAX = ONE/SAFMIN;
      ULP = SLAMCH( 'PRECISION' );
      SMLNUM = SAFMIN*( REAL( N )/ULP );

      BNORM = CLANHS( 'F', IHI-ILO+1, B( ILO, ILO ), LDB, RWORK );
      BTOL = max( SAFMIN, ULP*BNORM );

      ISTART = ILO;
      ISTOP = IHI;
      MAXIT = 30*( IHI-ILO+1 );
      LD = 0;

      for (IITER = 1; IITER <= MAXIT; IITER++) {
         if ( IITER >= MAXIT ) {
            INFO = ISTOP+1;
            GOTO 80;
         }
         if ( ISTART+1 >= ISTOP ) {
            ISTOP = ISTART;
            EXIT;
         }

         // Check deflations at the end
         if ( ABS( A( ISTOP, ISTOP-1 ) ) <= max( SMLNUM, ULP*( ABS( A( ISTOP, ISTOP ) )+ABS( A( ISTOP-1, ISTOP-1 ) ) ) ) ) {
            A( ISTOP, ISTOP-1 ) = CZERO;
            ISTOP = ISTOP-1;
            LD = 0;
            ESHIFT = CZERO;
         }
         // Check deflations at the start
         if ( ABS( A( ISTART+1, ISTART ) ) <= max( SMLNUM, ULP*( ABS( A( ISTART, ISTART ) )+ABS( A( ISTART+1, ISTART+1 ) ) ) ) ) {
            A( ISTART+1, ISTART ) = CZERO;
            ISTART = ISTART+1;
            LD = 0;
            ESHIFT = CZERO;
         }

         if ( ISTART+1 >= ISTOP ) {
            EXIT;
         }

         // Check interior deflations
         ISTART2 = ISTART;
         DO K = ISTOP, ISTART+1, -1;
            if ( ABS( A( K, K-1 ) ) <= max( SMLNUM, ULP*( ABS( A( K, K ) )+ABS( A( K-1, K-1 ) ) ) ) ) {
               A( K, K-1 ) = CZERO;
               ISTART2 = K;
               EXIT;
            }
         }

         // Get range to apply rotations to
         if ( ILSCHUR ) {
            ISTARTM = 1;
            ISTOPM = N;
         } else {
            ISTARTM = ISTART2;
            ISTOPM = ISTOP;
         }

         // Check infinite eigenvalues, this is done without blocking so might
         // slow down the method when many infinite eigenvalues are present
         K = ISTOP;
         DO WHILE ( K >= ISTART2 );

            if ( ABS( B( K, K ) ) < BTOL ) {
               // A diagonal element of B is negligible, move it
               // to the top and deflate it

               DO K2 = K, ISTART2+1, -1;
                  clartg(B( K2-1, K2 ), B( K2-1, K2-1 ), C1, S1, TEMP );
                  B( K2-1, K2 ) = TEMP;
                  B( K2-1, K2-1 ) = CZERO;
                   crot(K2-2-ISTARTM+1, B( ISTARTM, K2 ), 1, B( ISTARTM, K2-1 ), 1, C1, S1 );
                  crot(min( K2+1, ISTOP )-ISTARTM+1, A( ISTARTM, K2 ), 1, A( ISTARTM, K2-1 ), 1, C1, S1 );
                  if ( ILZ ) {
                     crot(N, Z( 1, K2 ), 1, Z( 1, K2-1 ), 1, C1, S1 );
                  }

                  if ( K2 < ISTOP ) {
                     clartg(A( K2, K2-1 ), A( K2+1, K2-1 ), C1, S1, TEMP );
                     A( K2, K2-1 ) = TEMP;
                     A( K2+1, K2-1 ) = CZERO;
                      crot(ISTOPM-K2+1, A( K2, K2 ), LDA, A( K2+1, K2 ), LDA, C1, S1 );
                     crot(ISTOPM-K2+1, B( K2, K2 ), LDB, B( K2+1, K2 ), LDB, C1, S1 );
                     if ( ILQ ) {
                        crot(N, Q( 1, K2 ), 1, Q( 1, K2+1 ), 1, C1, CONJG( S1 ) );
                     }
                  }

               }

               if ( ISTART2 < ISTOP ) {
                  clartg(A( ISTART2, ISTART2 ), A( ISTART2+1, ISTART2 ), C1, S1, TEMP );
                  A( ISTART2, ISTART2 ) = TEMP;
                  A( ISTART2+1, ISTART2 ) = CZERO;
                   crot(ISTOPM-( ISTART2+1 )+1, A( ISTART2, ISTART2+1 ), LDA, A( ISTART2+1, ISTART2+1 ), LDA, C1, S1 );
                  crot(ISTOPM-( ISTART2+1 )+1, B( ISTART2, ISTART2+1 ), LDB, B( ISTART2+1, ISTART2+1 ), LDB, C1, S1 );
                  if ( ILQ ) {
                     crot(N, Q( 1, ISTART2 ), 1, Q( 1, ISTART2+1 ), 1, C1, CONJG( S1 ) );
                  }
               }

               ISTART2 = ISTART2+1;

            }
            K = K-1;
         }

         // istart2 now points to the top of the bottom right
         // unreduced Hessenberg block
         if ( ISTART2 >= ISTOP ) {
            ISTOP = ISTART2-1;
            LD = 0;
            ESHIFT = CZERO;
            CYCLE;
         }

         NW = NWR;
         NSHIFTS = NSR;
         NBLOCK = NBR;

         if ( ISTOP-ISTART2+1 < NMIN ) {
            // Setting nw to the size of the subblock will make AED deflate
            // all the eigenvalues. This is slightly more efficient than just
            // using CHGEQZ because the off diagonal part gets updated via BLAS.
            if ( ISTOP-ISTART+1 < NMIN ) {
               NW = ISTOP-ISTART+1;
               ISTART2 = ISTART;
            } else {
               NW = ISTOP-ISTART2+1;
            }
         }


         // Time for AED

         claqz2(ILSCHUR, ILQ, ILZ, N, ISTART2, ISTOP, NW, A, LDA, B, LDB, Q, LDQ, Z, LDZ, N_UNDEFLATED, N_DEFLATED, ALPHA, BETA, WORK, NW, WORK( NW**2+1 ), NW, WORK( 2*NW**2+1 ), LWORK-2*NW**2, RWORK, REC, AED_INFO );

         if ( N_DEFLATED > 0 ) {
            ISTOP = ISTOP-N_DEFLATED;
            LD = 0;
            ESHIFT = CZERO;
         }
          if ( 100*N_DEFLATED > NIBBLE*( N_DEFLATED+N_UNDEFLATED ) || ISTOP-ISTART2+1 < NMIN ) {
            // AED has uncovered many eigenvalues. Skip a QZ sweep and run
            // AED again.
            CYCLE;
         }

         LD = LD+1;

         NS = min( NSHIFTS, ISTOP-ISTART2 );
         NS = min( NS, N_UNDEFLATED );
         SHIFTPOS = ISTOP-N_UNDEFLATED+1;

         if ( MOD( LD, 6 ) == 0 ) {

            // Exceptional shift.  Chosen for no particularly good reason.

            if ( ( REAL( MAXIT )*SAFMIN )*ABS( A( ISTOP, ISTOP-1 ) ) < ABS( A( ISTOP-1, ISTOP-1 ) ) ) {
               ESHIFT = A( ISTOP, ISTOP-1 )/B( ISTOP-1, ISTOP-1 );
            } else {
               ESHIFT = ESHIFT+CONE/( SAFMIN*REAL( MAXIT ) );
            }
            ALPHA( SHIFTPOS ) = CONE;
            BETA( SHIFTPOS ) = ESHIFT;
            NS = 1;
         }


         // Time for a QZ sweep

         claqz3(ILSCHUR, ILQ, ILZ, N, ISTART2, ISTOP, NS, NBLOCK, ALPHA( SHIFTPOS ), BETA( SHIFTPOS ), A, LDA, B, LDB, Q, LDQ, Z, LDZ, WORK, NBLOCK, WORK( NBLOCK** 2+1 ), NBLOCK, WORK( 2*NBLOCK**2+1 ), LWORK-2*NBLOCK**2, SWEEP_INFO );

      }


      // Call CHGEQZ to normalize the eigenvalue blocks and set the eigenvalues
      // If all the eigenvalues have been found, CHGEQZ will not do any iterations
      // and only normalize the blocks. In case of a rare convergence failure,
      // the single shift might perform better.

   80 CALL CHGEQZ( WANTS, WANTQ, WANTZ, N, ILO, IHI, A, LDA, B, LDB, ALPHA, BETA, Q, LDQ, Z, LDZ, WORK, LWORK, RWORK, NORM_INFO );

      INFO = NORM_INFO;

      END SUBROUTINE;
