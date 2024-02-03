      RECURSIVE SUBROUTINE SLAQZ0( WANTS, WANTQ, WANTZ, N, ILO, IHI, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, Q, LDQ, Z, LDZ, WORK, LWORK, REC, INFO )
      IMPLICIT NONE

      // Arguments
      String   , INTENT( IN ) :: WANTS, WANTQ, WANTZ;
      int    , INTENT( IN ) :: N, ILO, IHI, LDA, LDB, LDQ, LDZ, LWORK, REC;

      int    , INTENT( OUT ) :: INFO;
       REAL, INTENT( INOUT ) :: A( LDA, * ), B( LDB, * ), Q( LDQ, * ), Z( LDZ, * ), ALPHAR( * ), ALPHAI( * ), BETA( * ), WORK( * )

      // Parameters
      REAL :: ZERO, ONE, HALF
      const    ZERO = 0.0, ONE = 1.0, HALF = 0.5 ;

      // Local scalars
      REAL :: SMLNUM, ULP, ESHIFT, SAFMIN, SAFMAX, C1, S1, TEMP, SWAP, BNORM, BTOL;
      int     :: ISTART, ISTOP, IITER, MAXIT, ISTART2, K, LD, NSHIFTS, NBLOCK, NW, NMIN, NIBBLE, N_UNDEFLATED, N_DEFLATED, NS, SWEEP_INFO, SHIFTPOS, LWORKREQ, K2, ISTARTM, ISTOPM, IWANTS, IWANTQ, IWANTZ, NORM_INFO, AED_INFO, NWR, NBR, NSR, ITEMP1, ITEMP2, RCOST, I;
      bool    :: ILSCHUR, ILQ, ILZ;
      String    :: JBCMPZ*3;

      // External Functions
      // EXTERNAL :: XERBLA, SHGEQZ, SLAQZ3, SLAQZ4, SLASET, SLARTG, SROT
      REAL, EXTERNAL :: SLAMCH, SLANHS, SROUNDUP_LWORK
      bool   , EXTERNAL :: LSAME;
      int    , EXTERNAL :: ILAENV;


      // Decode wantS,wantQ,wantZ

      if ( LSAME( WANTS, 'E' ) ) {
         ILSCHUR = .FALSE.
         IWANTS = 1
      } else if ( LSAME( WANTS, 'S' ) ) {
         ILSCHUR = .TRUE.
         IWANTS = 2
      } else {
         IWANTS = 0
      }

      if ( LSAME( WANTQ, 'N' ) ) {
         ILQ = .FALSE.
         IWANTQ = 1
      } else if ( LSAME( WANTQ, 'V' ) ) {
         ILQ = .TRUE.
         IWANTQ = 2
      } else if ( LSAME( WANTQ, 'I' ) ) {
         ILQ = .TRUE.
         IWANTQ = 3
      } else {
         IWANTQ = 0
      }

      if ( LSAME( WANTZ, 'N' ) ) {
         ILZ = .FALSE.
         IWANTZ = 1
      } else if ( LSAME( WANTZ, 'V' ) ) {
         ILZ = .TRUE.
         IWANTZ = 2
      } else if ( LSAME( WANTZ, 'I' ) ) {
         ILZ = .TRUE.
         IWANTZ = 3
      } else {
         IWANTZ = 0
      }

      // Check Argument Values

      INFO = 0
      if ( IWANTS.EQ.0 ) {
         INFO = -1
      } else if ( IWANTQ.EQ.0 ) {
         INFO = -2
      } else if ( IWANTZ.EQ.0 ) {
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
         xerbla('SLAQZ0', -INFO );
         RETURN
      }


      // Quick return if possible

      if ( N.LE.0 ) {
         WORK( 1 ) = REAL( 1 )
         RETURN
      }


      // Get the parameters

      JBCMPZ( 1:1 ) = WANTS
      JBCMPZ( 2:2 ) = WANTQ
      JBCMPZ( 3:3 ) = WANTZ

      NMIN = ILAENV( 12, 'SLAQZ0', JBCMPZ, N, ILO, IHI, LWORK )

      NWR = ILAENV( 13, 'SLAQZ0', JBCMPZ, N, ILO, IHI, LWORK )
      NWR = MAX( 2, NWR )
      NWR = MIN( IHI-ILO+1, ( N-1 ) / 3, NWR )

      NIBBLE = ILAENV( 14, 'SLAQZ0', JBCMPZ, N, ILO, IHI, LWORK )

      NSR = ILAENV( 15, 'SLAQZ0', JBCMPZ, N, ILO, IHI, LWORK )
      NSR = MIN( NSR, ( N+6 ) / 9, IHI-ILO )
      NSR = MAX( 2, NSR-MOD( NSR, 2 ) )

      RCOST = ILAENV( 17, 'SLAQZ0', JBCMPZ, N, ILO, IHI, LWORK )
      ITEMP1 = INT( NSR/SQRT( 1+2*NSR/( REAL( RCOST )/100*N ) ) )
      ITEMP1 = ( ( ITEMP1-1 )/4 )*4+4
      NBR = NSR+ITEMP1

      if ( N .LT. NMIN .OR. REC .GE. 2 ) {
         shgeqz(WANTS, WANTQ, WANTZ, N, ILO, IHI, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, Q, LDQ, Z, LDZ, WORK, LWORK, INFO );
         RETURN
      }


      // Find out required workspace


      // Workspace query to slaqz3
      NW = MAX( NWR, NMIN )
      slaqz3(ILSCHUR, ILQ, ILZ, N, ILO, IHI, NW, A, LDA, B, LDB, Q, LDQ, Z, LDZ, N_UNDEFLATED, N_DEFLATED, ALPHAR, ALPHAI, BETA, WORK, NW, WORK, NW, WORK, -1, REC, AED_INFO );
      ITEMP1 = INT( WORK( 1 ) )
      // Workspace query to slaqz4
      slaqz4(ILSCHUR, ILQ, ILZ, N, ILO, IHI, NSR, NBR, ALPHAR, ALPHAI, BETA, A, LDA, B, LDB, Q, LDQ, Z, LDZ, WORK, NBR, WORK, NBR, WORK, -1, SWEEP_INFO );
      ITEMP2 = INT( WORK( 1 ) )

      LWORKREQ = MAX( ITEMP1+2*NW**2, ITEMP2+2*NBR**2 )
      if ( LWORK .EQ.-1 ) {
         WORK( 1 ) = SROUNDUP_LWORK( LWORKREQ )
         RETURN
      } else if ( LWORK .LT. LWORKREQ ) {
         INFO = -19
      }
      if ( INFO.NE.0 ) {
         xerbla('SLAQZ0', INFO );
         RETURN
      }

      // Initialize Q and Z

      IF( IWANTQ.EQ.3 ) CALL SLASET( 'FULL', N, N, ZERO, ONE, Q, LDQ )
      IF( IWANTZ.EQ.3 ) CALL SLASET( 'FULL', N, N, ZERO, ONE, Z, LDZ )

      // Get machine constants
      SAFMIN = SLAMCH( 'SAFE MINIMUM' )
      SAFMAX = ONE/SAFMIN
      ULP = SLAMCH( 'PRECISION' )
      SMLNUM = SAFMIN*( REAL( N )/ULP )

      BNORM = SLANHS( 'F', IHI-ILO+1, B( ILO, ILO ), LDB, WORK )
      BTOL = MAX( SAFMIN, ULP*BNORM )

      ISTART = ILO
      ISTOP = IHI
      MAXIT = 3*( IHI-ILO+1 )
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
         if ( ABS( A( ISTOP-1, ISTOP-2 ) ) .LE. MAX( SMLNUM, ULP*( ABS( A( ISTOP-1, ISTOP-1 ) )+ABS( A( ISTOP-2, ISTOP-2 ) ) ) ) ) {
            A( ISTOP-1, ISTOP-2 ) = ZERO
            ISTOP = ISTOP-2
            LD = 0
            ESHIFT = ZERO
         } else if ( ABS( A( ISTOP, ISTOP-1 ) ) .LE. MAX( SMLNUM, ULP*( ABS( A( ISTOP, ISTOP ) )+ABS( A( ISTOP-1, ISTOP-1 ) ) ) ) ) {
            A( ISTOP, ISTOP-1 ) = ZERO
            ISTOP = ISTOP-1
            LD = 0
            ESHIFT = ZERO
         }
         // Check deflations at the start
         if ( ABS( A( ISTART+2, ISTART+1 ) ) .LE. MAX( SMLNUM, ULP*( ABS( A( ISTART+1, ISTART+1 ) )+ABS( A( ISTART+2, ISTART+2 ) ) ) ) ) {
            A( ISTART+2, ISTART+1 ) = ZERO
            ISTART = ISTART+2
            LD = 0
            ESHIFT = ZERO
         } else if ( ABS( A( ISTART+1, ISTART ) ) .LE. MAX( SMLNUM, ULP*( ABS( A( ISTART, ISTART ) )+ABS( A( ISTART+1, ISTART+1 ) ) ) ) ) {
            A( ISTART+1, ISTART ) = ZERO
            ISTART = ISTART+1
            LD = 0
            ESHIFT = ZERO
         }

         if ( ISTART+1 .GE. ISTOP ) {
            EXIT
         }

         // Check interior deflations
         ISTART2 = ISTART
         DO K = ISTOP, ISTART+1, -1
            if ( ABS( A( K, K-1 ) ) .LE. MAX( SMLNUM, ULP*( ABS( A( K, K ) )+ABS( A( K-1, K-1 ) ) ) ) ) {
               A( K, K-1 ) = ZERO
               ISTART2 = K
               EXIT
            }
         END DO

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
                  slartg(B( K2-1, K2 ), B( K2-1, K2-1 ), C1, S1, TEMP );
                  B( K2-1, K2 ) = TEMP
                  B( K2-1, K2-1 ) = ZERO
                   srot(K2-2-ISTARTM+1, B( ISTARTM, K2 ), 1, B( ISTARTM, K2-1 ), 1, C1, S1 )                   CALL SROT( MIN( K2+1, ISTOP )-ISTARTM+1, A( ISTARTM, K2 ), 1, A( ISTARTM, K2-1 ), 1, C1, S1 );
                  if ( ILZ ) {
                     srot(N, Z( 1, K2 ), 1, Z( 1, K2-1 ), 1, C1, S1 );
                  }

                  if ( K2.LT.ISTOP ) {
                     slartg(A( K2, K2-1 ), A( K2+1, K2-1 ), C1, S1, TEMP );
                     A( K2, K2-1 ) = TEMP
                     A( K2+1, K2-1 ) = ZERO
                      srot(ISTOPM-K2+1, A( K2, K2 ), LDA, A( K2+1, K2 ), LDA, C1, S1 )                      CALL SROT( ISTOPM-K2+1, B( K2, K2 ), LDB, B( K2+1, K2 ), LDB, C1, S1 );
                     if ( ILQ ) {
                        srot(N, Q( 1, K2 ), 1, Q( 1, K2+1 ), 1, C1, S1 );
                     }
                  }

               END DO

               if ( ISTART2.LT.ISTOP ) {
                  slartg(A( ISTART2, ISTART2 ), A( ISTART2+1, ISTART2 ), C1, S1, TEMP );
                  A( ISTART2, ISTART2 ) = TEMP
                  A( ISTART2+1, ISTART2 ) = ZERO
                   srot(ISTOPM-( ISTART2+1 )+1, A( ISTART2, ISTART2+1 ), LDA, A( ISTART2+1, ISTART2+1 ), LDA, C1, S1 )                   CALL SROT( ISTOPM-( ISTART2+1 )+1, B( ISTART2, ISTART2+1 ), LDB, B( ISTART2+1, ISTART2+1 ), LDB, C1, S1 );
                  if ( ILQ ) {
                     srot(N, Q( 1, ISTART2 ), 1, Q( 1, ISTART2+1 ), 1, C1, S1 );
                  }
               }

               ISTART2 = ISTART2+1

            }
            K = K-1
         END DO

         // istart2 now points to the top of the bottom right
         // unreduced Hessenberg block
         if ( ISTART2 .GE. ISTOP ) {
            ISTOP = ISTART2-1
            LD = 0
            ESHIFT = ZERO
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

         slaqz3(ILSCHUR, ILQ, ILZ, N, ISTART2, ISTOP, NW, A, LDA, B, LDB, Q, LDQ, Z, LDZ, N_UNDEFLATED, N_DEFLATED, ALPHAR, ALPHAI, BETA, WORK, NW, WORK( NW**2+1 ), NW, WORK( 2*NW**2+1 ), LWORK-2*NW**2, REC, AED_INFO );

         if ( N_DEFLATED > 0 ) {
            ISTOP = ISTOP-N_DEFLATED
            LD = 0
            ESHIFT = ZERO
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

         // Shuffle shifts to put double shifts in front
         // This ensures that we don't split up a double shift

         DO I = SHIFTPOS, SHIFTPOS+N_UNDEFLATED-1, 2
            if ( ALPHAI( I ).NE.-ALPHAI( I+1 ) ) {

               SWAP = ALPHAR( I )
               ALPHAR( I ) = ALPHAR( I+1 )
               ALPHAR( I+1 ) = ALPHAR( I+2 )
               ALPHAR( I+2 ) = SWAP

               SWAP = ALPHAI( I )
               ALPHAI( I ) = ALPHAI( I+1 )
               ALPHAI( I+1 ) = ALPHAI( I+2 )
               ALPHAI( I+2 ) = SWAP

               SWAP = BETA( I )
               BETA( I ) = BETA( I+1 )
               BETA( I+1 ) = BETA( I+2 )
               BETA( I+2 ) = SWAP
            }
         END DO

         if ( MOD( LD, 6 ) .EQ. 0 ) {

            // Exceptional shift.  Chosen for no particularly good reason.

            if ( ( REAL( MAXIT )*SAFMIN )*ABS( A( ISTOP, ISTOP-1 ) ).LT.ABS( A( ISTOP-1, ISTOP-1 ) ) ) {
               ESHIFT = A( ISTOP, ISTOP-1 )/B( ISTOP-1, ISTOP-1 )
            } else {
               ESHIFT = ESHIFT+ONE/( SAFMIN*REAL( MAXIT ) )
            }
            ALPHAR( SHIFTPOS ) = ONE
            ALPHAR( SHIFTPOS+1 ) = ZERO
            ALPHAI( SHIFTPOS ) = ZERO
            ALPHAI( SHIFTPOS+1 ) = ZERO
            BETA( SHIFTPOS ) = ESHIFT
            BETA( SHIFTPOS+1 ) = ESHIFT
            NS = 2
         }


         // Time for a QZ sweep

         slaqz4(ILSCHUR, ILQ, ILZ, N, ISTART2, ISTOP, NS, NBLOCK, ALPHAR( SHIFTPOS ), ALPHAI( SHIFTPOS ), BETA( SHIFTPOS ), A, LDA, B, LDB, Q, LDQ, Z, LDZ, WORK, NBLOCK, WORK( NBLOCK**2+1 ), NBLOCK, WORK( 2*NBLOCK**2+1 ), LWORK-2*NBLOCK**2, SWEEP_INFO );

      END DO


      // Call SHGEQZ to normalize the eigenvalue blocks and set the eigenvalues
      // If all the eigenvalues have been found, SHGEQZ will not do any iterations
      // and only normalize the blocks. In case of a rare convergence failure,
      // the single shift might perform better.

   80 CALL SHGEQZ( WANTS, WANTQ, WANTZ, N, ILO, IHI, A, LDA, B, LDB, ALPHAR, ALPHAI, BETA, Q, LDQ, Z, LDZ, WORK, LWORK, NORM_INFO )

      INFO = NORM_INFO

      END SUBROUTINE
