      RECURSIVE SUBROUTINE ZLAQZ0( WANTS, WANTQ, WANTZ, N, ILO, IHI, A, LDA, B, LDB, ALPHA, BETA, Q, LDQ, Z, LDZ, WORK, LWORK, RWORK, REC, INFO )
      IMPLICIT NONE

      // Arguments
      String   , INTENT( IN ) :: WANTS, WANTQ, WANTZ;
      int    , INTENT( IN ) :: N, ILO, IHI, LDA, LDB, LDQ, LDZ, LWORK, REC;
      int    , INTENT( OUT ) :: INFO;
      COMPLEX*16, INTENT( INOUT ) :: A( LDA, * ), B( LDB, * ), Q( LDQ, * ), Z( LDZ, * ), ALPHA( * ), BETA( * ), WORK( * )
      double          , INTENT( OUT ) :: RWORK( * );

      // Parameters
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) ;
      double           :: ZERO, ONE, HALF;
      const    ZERO = 0.0D0, ONE = 1.0D0, HALF = 0.5D0 ;

      // Local scalars
      double           :: SMLNUM, ULP, SAFMIN, SAFMAX, C1, TEMPR, BNORM, BTOL;
      COMPLEX*16 :: ESHIFT, S1, TEMP
      int     :: ISTART, ISTOP, IITER, MAXIT, ISTART2, K, LD, NSHIFTS, NBLOCK, NW, NMIN, NIBBLE, N_UNDEFLATED, N_DEFLATED, NS, SWEEP_INFO, SHIFTPOS, LWORKREQ, K2, ISTARTM, ISTOPM, IWANTS, IWANTQ, IWANTZ, NORM_INFO, AED_INFO, NWR, NBR, NSR, ITEMP1, ITEMP2, RCOST;
      bool    :: ILSCHUR, ILQ, ILZ;
      String    :: JBCMPZ*3;

      // External Functions
      // EXTERNAL :: XERBLA, ZHGEQZ, ZLAQZ2, ZLAQZ3, ZLASET, ZLARTG, ZROT
      double          , EXTERNAL :: DLAMCH, ZLANHS;
      bool   , EXTERNAL :: LSAME;
      int    , EXTERNAL :: ILAENV;


      // Decode wantS,wantQ,wantZ

      if ( LSAME( WANTS, 'E' ) ) {
         ILSCHUR = false;
         IWANTS = 1
      } else if ( LSAME( WANTS, 'S' ) ) {
         ILSCHUR = true;
         IWANTS = 2
      } else {
         IWANTS = 0
      }

      if ( LSAME( WANTQ, 'N' ) ) {
         ILQ = false;
         IWANTQ = 1
      } else if ( LSAME( WANTQ, 'V' ) ) {
         ILQ = true;
         IWANTQ = 2
      } else if ( LSAME( WANTQ, 'I' ) ) {
         ILQ = true;
         IWANTQ = 3
      } else {
         IWANTQ = 0
      }

      if ( LSAME( WANTZ, 'N' ) ) {
         ILZ = false;
         IWANTZ = 1
      } else if ( LSAME( WANTZ, 'V' ) ) {
         ILZ = true;
         IWANTZ = 2
      } else if ( LSAME( WANTZ, 'I' ) ) {
         ILZ = true;
         IWANTZ = 3
      } else {
         IWANTZ = 0
      }

      // Check Argument Values

      INFO = 0
      if ( IWANTS == 0 ) {
         INFO = -1
      } else if ( IWANTQ == 0 ) {
         INFO = -2
      } else if ( IWANTZ == 0 ) {
         INFO = -3
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( ILO.LT.1 ) {
         INFO = -5
      } else if ( IHI.GT.N .OR. IHI.LT.ILO-1 ) {
         INFO = -6
      } else if ( LDA.LT.N ) {
         INFO = -8
      } else if ( LDB.LT.N ) {
         INFO = -10
      } else if ( LDQ.LT.1 .OR. ( ILQ .AND. LDQ.LT.N ) ) {
         INFO = -15
      } else if ( LDZ.LT.1 .OR. ( ILZ .AND. LDZ.LT.N ) ) {
         INFO = -17
      }
      if ( INFO.NE.0 ) {
         xerbla('ZLAQZ0', -INFO );
         RETURN
      }


      // Quick return if possible

      if ( N.LE.0 ) {
         WORK( 1 ) = DBLE( 1 )
         RETURN
      }


      // Get the parameters

      JBCMPZ( 1:1 ) = WANTS
      JBCMPZ( 2:2 ) = WANTQ
      JBCMPZ( 3:3 ) = WANTZ

      NMIN = ILAENV( 12, 'ZLAQZ0', JBCMPZ, N, ILO, IHI, LWORK )

      NWR = ILAENV( 13, 'ZLAQZ0', JBCMPZ, N, ILO, IHI, LWORK )
      NWR = MAX( 2, NWR )
      NWR = MIN( IHI-ILO+1, ( N-1 ) / 3, NWR )

      NIBBLE = ILAENV( 14, 'ZLAQZ0', JBCMPZ, N, ILO, IHI, LWORK )

      NSR = ILAENV( 15, 'ZLAQZ0', JBCMPZ, N, ILO, IHI, LWORK )
      NSR = MIN( NSR, ( N+6 ) / 9, IHI-ILO )
      NSR = MAX( 2, NSR-MOD( NSR, 2 ) )

      RCOST = ILAENV( 17, 'ZLAQZ0', JBCMPZ, N, ILO, IHI, LWORK )
      ITEMP1 = INT( NSR/SQRT( 1+2*NSR/( DBLE( RCOST )/100*N ) ) )
      ITEMP1 = ( ( ITEMP1-1 )/4 )*4+4
      NBR = NSR+ITEMP1

      if ( N .LT. NMIN .OR. REC .GE. 2 ) {
         zhgeqz(WANTS, WANTQ, WANTZ, N, ILO, IHI, A, LDA, B, LDB, ALPHA, BETA, Q, LDQ, Z, LDZ, WORK, LWORK, RWORK, INFO );
         RETURN
      }


      // Find out required workspace


      // Workspace query to ZLAQZ2
      NW = MAX( NWR, NMIN )
      zlaqz2(ILSCHUR, ILQ, ILZ, N, ILO, IHI, NW, A, LDA, B, LDB, Q, LDQ, Z, LDZ, N_UNDEFLATED, N_DEFLATED, ALPHA, BETA, WORK, NW, WORK, NW, WORK, -1, RWORK, REC, AED_INFO );
      ITEMP1 = INT( WORK( 1 ) )
      // Workspace query to ZLAQZ3
      zlaqz3(ILSCHUR, ILQ, ILZ, N, ILO, IHI, NSR, NBR, ALPHA, BETA, A, LDA, B, LDB, Q, LDQ, Z, LDZ, WORK, NBR, WORK, NBR, WORK, -1, SWEEP_INFO );
      ITEMP2 = INT( WORK( 1 ) )

      LWORKREQ = MAX( ITEMP1+2*NW**2, ITEMP2+2*NBR**2 )
      if ( LWORK == -1 ) {
         WORK( 1 ) = DBLE( LWORKREQ )
         RETURN
      } else if ( LWORK .LT. LWORKREQ ) {
         INFO = -19
      }
      if ( INFO.NE.0 ) {
         xerbla('ZLAQZ0', INFO );
         RETURN
      }

      // Initialize Q and Z

      if (IWANTQ == 3) CALL ZLASET( 'FULL', N, N, CZERO, CONE, Q, LDQ )       IF( IWANTZ == 3 ) CALL ZLASET( 'FULL', N, N, CZERO, CONE, Z, LDZ );

      // Get machine constants
      SAFMIN = DLAMCH( 'SAFE MINIMUM' )
      SAFMAX = ONE/SAFMIN
      ULP = DLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( DBLE( N )/ULP )

      BNORM = ZLANHS( 'F', IHI-ILO+1, B( ILO, ILO ), LDB, RWORK )
      BTOL = MAX( SAFMIN, ULP*BNORM )

      ISTART = ILO
      ISTOP = IHI
      MAXIT = 30*( IHI-ILO+1 )
      LD = 0

      for (IITER = 1; IITER <= MAXIT; IITER++) {
         if ( IITER .GE. MAXIT ) {
            INFO = ISTOP+1
            GOTO 80
         }
         if ( ISTART+1 .GE. ISTOP ) {
            ISTOP = ISTART
            EXIT
         }

         // Check deflations at the end
         if ( ABS( A( ISTOP, ISTOP-1 ) ) .LE. MAX( SMLNUM, ULP*( ABS( A( ISTOP, ISTOP ) )+ABS( A( ISTOP-1, ISTOP-1 ) ) ) ) ) {
            A( ISTOP, ISTOP-1 ) = CZERO
            ISTOP = ISTOP-1
            LD = 0
            ESHIFT = CZERO
         }
         // Check deflations at the start
         if ( ABS( A( ISTART+1, ISTART ) ) .LE. MAX( SMLNUM, ULP*( ABS( A( ISTART, ISTART ) )+ABS( A( ISTART+1, ISTART+1 ) ) ) ) ) {
            A( ISTART+1, ISTART ) = CZERO
            ISTART = ISTART+1
            LD = 0
            ESHIFT = CZERO
         }

         if ( ISTART+1 .GE. ISTOP ) {
            EXIT
         }

         // Check interior deflations
         ISTART2 = ISTART
         DO K = ISTOP, ISTART+1, -1
            if ( ABS( A( K, K-1 ) ) .LE. MAX( SMLNUM, ULP*( ABS( A( K, K ) )+ABS( A( K-1, K-1 ) ) ) ) ) {
               A( K, K-1 ) = CZERO
               ISTART2 = K
               EXIT
            }
         }

         // Get range to apply rotations to
         if ( ILSCHUR ) {
            ISTARTM = 1
            ISTOPM = N
         } else {
            ISTARTM = ISTART2
            ISTOPM = ISTOP
         }

         // Check infinite eigenvalues, this is done without blocking so might
         // slow down the method when many infinite eigenvalues are present
         K = ISTOP
         DO WHILE ( K.GE.ISTART2 )

            if ( ABS( B( K, K ) ) .LT. BTOL ) {
               // A diagonal element of B is negligible, move it
               // to the top and deflate it

               DO K2 = K, ISTART2+1, -1
                  zlartg(B( K2-1, K2 ), B( K2-1, K2-1 ), C1, S1, TEMP );
                  B( K2-1, K2 ) = TEMP
                  B( K2-1, K2-1 ) = CZERO
                   zrot(K2-2-ISTARTM+1, B( ISTARTM, K2 ), 1, B( ISTARTM, K2-1 ), 1, C1, S1 );
                  zrot(MIN( K2+1, ISTOP )-ISTARTM+1, A( ISTARTM, K2 ), 1, A( ISTARTM, K2-1 ), 1, C1, S1 );
                  if ( ILZ ) {
                     zrot(N, Z( 1, K2 ), 1, Z( 1, K2-1 ), 1, C1, S1 );
                  }

                  if ( K2.LT.ISTOP ) {
                     zlartg(A( K2, K2-1 ), A( K2+1, K2-1 ), C1, S1, TEMP );
                     A( K2, K2-1 ) = TEMP
                     A( K2+1, K2-1 ) = CZERO
                      zrot(ISTOPM-K2+1, A( K2, K2 ), LDA, A( K2+1, K2 ), LDA, C1, S1 );
                     zrot(ISTOPM-K2+1, B( K2, K2 ), LDB, B( K2+1, K2 ), LDB, C1, S1 );
                     if ( ILQ ) {
                        zrot(N, Q( 1, K2 ), 1, Q( 1, K2+1 ), 1, C1, DCONJG( S1 ) );
                     }
                  }

               }

               if ( ISTART2.LT.ISTOP ) {
                  zlartg(A( ISTART2, ISTART2 ), A( ISTART2+1, ISTART2 ), C1, S1, TEMP );
                  A( ISTART2, ISTART2 ) = TEMP
                  A( ISTART2+1, ISTART2 ) = CZERO
                   zrot(ISTOPM-( ISTART2+1 )+1, A( ISTART2, ISTART2+1 ), LDA, A( ISTART2+1, ISTART2+1 ), LDA, C1, S1 );
                  zrot(ISTOPM-( ISTART2+1 )+1, B( ISTART2, ISTART2+1 ), LDB, B( ISTART2+1, ISTART2+1 ), LDB, C1, S1 );
                  if ( ILQ ) {
                     zrot(N, Q( 1, ISTART2 ), 1, Q( 1, ISTART2+1 ), 1, C1, DCONJG( S1 ) );
                  }
               }

               ISTART2 = ISTART2+1

            }
            K = K-1
         }

         // istart2 now points to the top of the bottom right
         // unreduced Hessenberg block
         if ( ISTART2 .GE. ISTOP ) {
            ISTOP = ISTART2-1
            LD = 0
            ESHIFT = CZERO
            CYCLE
         }

         NW = NWR
         NSHIFTS = NSR
         NBLOCK = NBR

         if ( ISTOP-ISTART2+1 .LT. NMIN ) {
            // Setting nw to the size of the subblock will make AED deflate
            // all the eigenvalues. This is slightly more efficient than just
            // using qz_small because the off diagonal part gets updated via BLAS.
            if ( ISTOP-ISTART+1 .LT. NMIN ) {
               NW = ISTOP-ISTART+1
               ISTART2 = ISTART
            } else {
               NW = ISTOP-ISTART2+1
            }
         }


         // Time for AED

         zlaqz2(ILSCHUR, ILQ, ILZ, N, ISTART2, ISTOP, NW, A, LDA, B, LDB, Q, LDQ, Z, LDZ, N_UNDEFLATED, N_DEFLATED, ALPHA, BETA, WORK, NW, WORK( NW**2+1 ), NW, WORK( 2*NW**2+1 ), LWORK-2*NW**2, RWORK, REC, AED_INFO );

         if ( N_DEFLATED > 0 ) {
            ISTOP = ISTOP-N_DEFLATED
            LD = 0
            ESHIFT = CZERO
         }
          if ( 100*N_DEFLATED > NIBBLE*( N_DEFLATED+N_UNDEFLATED ) .OR. ISTOP-ISTART2+1 .LT. NMIN ) {
            // AED has uncovered many eigenvalues. Skip a QZ sweep and run
            // AED again.
            CYCLE
         }

         LD = LD+1

         NS = MIN( NSHIFTS, ISTOP-ISTART2 )
         NS = MIN( NS, N_UNDEFLATED )
         SHIFTPOS = ISTOP-N_UNDEFLATED+1

         if ( MOD( LD, 6 ) == 0 ) {

            // Exceptional shift.  Chosen for no particularly good reason.

            if ( ( DBLE( MAXIT )*SAFMIN )*ABS( A( ISTOP, ISTOP-1 ) ).LT.ABS( A( ISTOP-1, ISTOP-1 ) ) ) {
               ESHIFT = A( ISTOP, ISTOP-1 )/B( ISTOP-1, ISTOP-1 )
            } else {
               ESHIFT = ESHIFT+CONE/( SAFMIN*DBLE( MAXIT ) )
            }
            ALPHA( SHIFTPOS ) = CONE
            BETA( SHIFTPOS ) = ESHIFT
            NS = 1
         }


         // Time for a QZ sweep

         zlaqz3(ILSCHUR, ILQ, ILZ, N, ISTART2, ISTOP, NS, NBLOCK, ALPHA( SHIFTPOS ), BETA( SHIFTPOS ), A, LDA, B, LDB, Q, LDQ, Z, LDZ, WORK, NBLOCK, WORK( NBLOCK** 2+1 ), NBLOCK, WORK( 2*NBLOCK**2+1 ), LWORK-2*NBLOCK**2, SWEEP_INFO );

      }


      // Call ZHGEQZ to normalize the eigenvalue blocks and set the eigenvalues
      // If all the eigenvalues have been found, ZHGEQZ will not do any iterations
      // and only normalize the blocks. In case of a rare convergence failure,
      // the single shift might perform better.

   80 CALL ZHGEQZ( WANTS, WANTQ, WANTZ, N, ILO, IHI, A, LDA, B, LDB, ALPHA, BETA, Q, LDQ, Z, LDZ, WORK, LWORK, RWORK, NORM_INFO )

      INFO = NORM_INFO

      END SUBROUTINE
