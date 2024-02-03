      SUBROUTINE SBDSDC( UPLO, COMPQ, N, D, E, U, LDU, VT, LDVT, Q, IQ, WORK, IWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             COMPQ, UPLO;
      int                INFO, LDU, LDVT, N;
      // ..
      // .. Array Arguments ..
      int                IQ( * ), IWORK( * );
      REAL               D( * ), E( * ), Q( * ), U( LDU, * ), VT( LDVT, * ), WORK( * )
      // ..

*  =====================================================================
*  Changed dimension statement in comment describing E from (N) to
*  (N-1).  Sven, 17 Feb 05.
*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TWO
      const              ZERO = 0.0E+0, ONE = 1.0E+0, TWO = 2.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                DIFL, DIFR, GIVCOL, GIVNUM, GIVPTR, I, IC, ICOMPQ, IERR, II, IS, IU, IUPLO, IVT, J, K, KK, MLVL, NM1, NSIZE, PERM, POLES, QSTART, SMLSIZ, SMLSZP, SQRE, START, WSTART, Z;
      REAL               CS, EPS, ORGNRM, P, R, SN
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               SLAMCH, SLANST
      // EXTERNAL SLAMCH, SLANST, ILAENV, LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SLARTG, SLASCL, SLASD0, SLASDA, SLASDQ, SLASET, SLASR, SSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL, ABS, INT, LOG, SIGN
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0

      IUPLO = 0
      IF( LSAME( UPLO, 'U' ) ) IUPLO = 1       IF( LSAME( UPLO, 'L' ) ) IUPLO = 2
      if ( LSAME( COMPQ, 'N' ) ) {
         ICOMPQ = 0
      } else if ( LSAME( COMPQ, 'P' ) ) {
         ICOMPQ = 1
      } else if ( LSAME( COMPQ, 'I' ) ) {
         ICOMPQ = 2
      } else {
         ICOMPQ = -1
      }
      if ( IUPLO.EQ.0 ) {
         INFO = -1
      } else if ( ICOMPQ.LT.0 ) {
         INFO = -2
      } else if ( N.LT.0 ) {
         INFO = -3
      } else if ( ( LDU.LT.1 ) .OR. ( ( ICOMPQ.EQ.2 ) .AND. ( LDU.LT. N ) ) ) {
         INFO = -7
      } else if ( ( LDVT.LT.1 ) .OR. ( ( ICOMPQ.EQ.2 ) .AND. ( LDVT.LT. N ) ) ) {
         INFO = -9
      }
      if ( INFO.NE.0 ) {
         xerbla('SBDSDC', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN
      SMLSIZ = ILAENV( 9, 'SBDSDC', ' ', 0, 0, 0, 0 )
      if ( N.EQ.1 ) {
         if ( ICOMPQ.EQ.1 ) {
            Q( 1 ) = SIGN( ONE, D( 1 ) )
            Q( 1+SMLSIZ*N ) = ONE
         } else if ( ICOMPQ.EQ.2 ) {
            U( 1, 1 ) = SIGN( ONE, D( 1 ) )
            VT( 1, 1 ) = ONE
         }
         D( 1 ) = ABS( D( 1 ) )
         RETURN
      }
      NM1 = N - 1

      // If matrix lower bidiagonal, rotate to be upper bidiagonal
      // by applying Givens rotations on the left

      WSTART = 1
      QSTART = 3
      if ( ICOMPQ.EQ.1 ) {
         scopy(N, D, 1, Q( 1 ), 1 );
         scopy(N-1, E, 1, Q( N+1 ), 1 );
      }
      if ( IUPLO.EQ.2 ) {
         QSTART = 5
         IF( ICOMPQ .EQ. 2 ) WSTART = 2*N - 1
         DO 10 I = 1, N - 1
            slartg(D( I ), E( I ), CS, SN, R );
            D( I ) = R
            E( I ) = SN*D( I+1 )
            D( I+1 ) = CS*D( I+1 )
            if ( ICOMPQ.EQ.1 ) {
               Q( I+2*N ) = CS
               Q( I+3*N ) = SN
            } else if ( ICOMPQ.EQ.2 ) {
               WORK( I ) = CS
               WORK( NM1+I ) = -SN
            }
         } // 10
      }

      // If ICOMPQ = 0, use SLASDQ to compute the singular values.

      if ( ICOMPQ.EQ.0 ) {
         // Ignore WSTART, instead using WORK( 1 ), since the two vectors
         // for CS and -SN above are added only if ICOMPQ == 2,
         // and adding them exceeds documented WORK size of 4*n.
         slasdq('U', 0, N, 0, 0, 0, D, E, VT, LDVT, U, LDU, U, LDU, WORK( 1 ), INFO );
         GO TO 40
      }

      // If N is smaller than the minimum divide size SMLSIZ, then solve
      // the problem with another solver.

      if ( N.LE.SMLSIZ ) {
         if ( ICOMPQ.EQ.2 ) {
            slaset('A', N, N, ZERO, ONE, U, LDU );
            slaset('A', N, N, ZERO, ONE, VT, LDVT );
            slasdq('U', 0, N, N, N, 0, D, E, VT, LDVT, U, LDU, U, LDU, WORK( WSTART ), INFO );
         } else if ( ICOMPQ.EQ.1 ) {
            IU = 1
            IVT = IU + N
            slaset('A', N, N, ZERO, ONE, Q( IU+( QSTART-1 )*N ), N )             CALL SLASET( 'A', N, N, ZERO, ONE, Q( IVT+( QSTART-1 )*N ), N )             CALL SLASDQ( 'U', 0, N, N, N, 0, D, E, Q( IVT+( QSTART-1 )*N ), N, Q( IU+( QSTART-1 )*N ), N, Q( IU+( QSTART-1 )*N ), N, WORK( WSTART ), INFO );
         }
         GO TO 40
      }

      if ( ICOMPQ.EQ.2 ) {
         slaset('A', N, N, ZERO, ONE, U, LDU );
         slaset('A', N, N, ZERO, ONE, VT, LDVT );
      }

      // Scale.

      ORGNRM = SLANST( 'M', N, D, E )
      IF( ORGNRM.EQ.ZERO ) RETURN
      slascl('G', 0, 0, ORGNRM, ONE, N, 1, D, N, IERR );
      slascl('G', 0, 0, ORGNRM, ONE, NM1, 1, E, NM1, IERR );

      EPS = SLAMCH( 'Epsilon' )

      MLVL = INT( LOG( REAL( N ) / REAL( SMLSIZ+1 ) ) / LOG( TWO ) ) + 1
      SMLSZP = SMLSIZ + 1

      if ( ICOMPQ.EQ.1 ) {
         IU = 1
         IVT = 1 + SMLSIZ
         DIFL = IVT + SMLSZP
         DIFR = DIFL + MLVL
         Z = DIFR + MLVL*2
         IC = Z + MLVL
         IS = IC + 1
         POLES = IS + 1
         GIVNUM = POLES + 2*MLVL

         K = 1
         GIVPTR = 2
         PERM = 3
         GIVCOL = PERM + MLVL
      }

      for (I = 1; I <= N; I++) { // 20
         if ( ABS( D( I ) ).LT.EPS ) {
            D( I ) = SIGN( EPS, D( I ) )
         }
      } // 20

      START = 1
      SQRE = 0

      for (I = 1; I <= NM1; I++) { // 30
         if ( ( ABS( E( I ) ).LT.EPS ) .OR. ( I.EQ.NM1 ) ) {

         // Subproblem found. First determine its size and then
         // apply divide and conquer on it.

            if ( I.LT.NM1 ) {

         // A subproblem with E(I) small for I < NM1.

               NSIZE = I - START + 1
            } else if ( ABS( E( I ) ).GE.EPS ) {

         // A subproblem with E(NM1) not too small but I = NM1.

               NSIZE = N - START + 1
            } else {

         // A subproblem with E(NM1) small. This implies an
         // 1-by-1 subproblem at D(N). Solve this 1-by-1 problem
         // first.

               NSIZE = I - START + 1
               if ( ICOMPQ.EQ.2 ) {
                  U( N, N ) = SIGN( ONE, D( N ) )
                  VT( N, N ) = ONE
               } else if ( ICOMPQ.EQ.1 ) {
                  Q( N+( QSTART-1 )*N ) = SIGN( ONE, D( N ) )
                  Q( N+( SMLSIZ+QSTART-1 )*N ) = ONE
               }
               D( N ) = ABS( D( N ) )
            }
            if ( ICOMPQ.EQ.2 ) {
               slasd0(NSIZE, SQRE, D( START ), E( START ), U( START, START ), LDU, VT( START, START ), LDVT, SMLSIZ, IWORK, WORK( WSTART ), INFO );
            } else {
               slasda(ICOMPQ, SMLSIZ, NSIZE, SQRE, D( START ), E( START ), Q( START+( IU+QSTART-2 )*N ), N, Q( START+( IVT+QSTART-2 )*N ), IQ( START+K*N ), Q( START+( DIFL+QSTART-2 )* N ), Q( START+( DIFR+QSTART-2 )*N ), Q( START+( Z+QSTART-2 )*N ), Q( START+( POLES+QSTART-2 )*N ), IQ( START+GIVPTR*N ), IQ( START+GIVCOL*N ), N, IQ( START+PERM*N ), Q( START+( GIVNUM+QSTART-2 )*N ), Q( START+( IC+QSTART-2 )*N ), Q( START+( IS+QSTART-2 )*N ), WORK( WSTART ), IWORK, INFO );
            }
            if ( INFO.NE.0 ) {
               RETURN
            }
            START = I + 1
         }
      } // 30

      // Unscale

      slascl('G', 0, 0, ONE, ORGNRM, N, 1, D, N, IERR );
      } // 40

      // Use Selection Sort to minimize swaps of singular vectors

      for (II = 2; II <= N; II++) { // 60
         I = II - 1
         KK = I
         P = D( I )
         for (J = II; J <= N; J++) { // 50
            if ( D( J ).GT.P ) {
               KK = J
               P = D( J )
            }
         } // 50
         if ( KK.NE.I ) {
            D( KK ) = D( I )
            D( I ) = P
            if ( ICOMPQ.EQ.1 ) {
               IQ( I ) = KK
            } else if ( ICOMPQ.EQ.2 ) {
               sswap(N, U( 1, I ), 1, U( 1, KK ), 1 );
               sswap(N, VT( I, 1 ), LDVT, VT( KK, 1 ), LDVT );
            }
         } else if ( ICOMPQ.EQ.1 ) {
            IQ( I ) = I
         }
      } // 60

      // If ICOMPQ = 1, use IQ(N,1) as the indicator for UPLO

      if ( ICOMPQ.EQ.1 ) {
         if ( IUPLO.EQ.1 ) {
            IQ( N ) = 1
         } else {
            IQ( N ) = 0
         }
      }

      // If B is lower bidiagonal, update U by those Givens rotations
      // which rotated B to be upper bidiagonal

      IF( ( IUPLO.EQ.2 ) .AND. ( ICOMPQ.EQ.2 ) ) CALL SLASR( 'L', 'V', 'B', N, N, WORK( 1 ), WORK( N ), U, LDU )

      RETURN

      // End of SBDSDC

      }
