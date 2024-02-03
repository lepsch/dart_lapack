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
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0, TWO = 2.0E+0 )
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
      IF( LSAME( COMPQ, 'N' ) ) THEN
         ICOMPQ = 0
      ELSE IF( LSAME( COMPQ, 'P' ) ) THEN
         ICOMPQ = 1
      ELSE IF( LSAME( COMPQ, 'I' ) ) THEN
         ICOMPQ = 2
      ELSE
         ICOMPQ = -1
      END IF
      IF( IUPLO.EQ.0 ) THEN
         INFO = -1
      ELSE IF( ICOMPQ.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( ( LDU.LT.1 ) .OR. ( ( ICOMPQ.EQ.2 ) .AND. ( LDU.LT. N ) ) ) THEN
         INFO = -7
      ELSE IF( ( LDVT.LT.1 ) .OR. ( ( ICOMPQ.EQ.2 ) .AND. ( LDVT.LT. N ) ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SBDSDC', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) RETURN
      SMLSIZ = ILAENV( 9, 'SBDSDC', ' ', 0, 0, 0, 0 )
      IF( N.EQ.1 ) THEN
         IF( ICOMPQ.EQ.1 ) THEN
            Q( 1 ) = SIGN( ONE, D( 1 ) )
            Q( 1+SMLSIZ*N ) = ONE
         ELSE IF( ICOMPQ.EQ.2 ) THEN
            U( 1, 1 ) = SIGN( ONE, D( 1 ) )
            VT( 1, 1 ) = ONE
         END IF
         D( 1 ) = ABS( D( 1 ) )
         RETURN
      END IF
      NM1 = N - 1

      // If matrix lower bidiagonal, rotate to be upper bidiagonal
      // by applying Givens rotations on the left

      WSTART = 1
      QSTART = 3
      IF( ICOMPQ.EQ.1 ) THEN
         CALL SCOPY( N, D, 1, Q( 1 ), 1 )
         CALL SCOPY( N-1, E, 1, Q( N+1 ), 1 )
      END IF
      IF( IUPLO.EQ.2 ) THEN
         QSTART = 5
         IF( ICOMPQ .EQ. 2 ) WSTART = 2*N - 1
         DO 10 I = 1, N - 1
            CALL SLARTG( D( I ), E( I ), CS, SN, R )
            D( I ) = R
            E( I ) = SN*D( I+1 )
            D( I+1 ) = CS*D( I+1 )
            IF( ICOMPQ.EQ.1 ) THEN
               Q( I+2*N ) = CS
               Q( I+3*N ) = SN
            ELSE IF( ICOMPQ.EQ.2 ) THEN
               WORK( I ) = CS
               WORK( NM1+I ) = -SN
            END IF
   10    CONTINUE
      END IF

      // If ICOMPQ = 0, use SLASDQ to compute the singular values.

      IF( ICOMPQ.EQ.0 ) THEN
         // Ignore WSTART, instead using WORK( 1 ), since the two vectors
         // for CS and -SN above are added only if ICOMPQ == 2,
         // and adding them exceeds documented WORK size of 4*n.
         CALL SLASDQ( 'U', 0, N, 0, 0, 0, D, E, VT, LDVT, U, LDU, U, LDU, WORK( 1 ), INFO )
         GO TO 40
      END IF

      // If N is smaller than the minimum divide size SMLSIZ, then solve
     t // he problem with another solver.

      IF( N.LE.SMLSIZ ) THEN
         IF( ICOMPQ.EQ.2 ) THEN
            CALL SLASET( 'A', N, N, ZERO, ONE, U, LDU )
            CALL SLASET( 'A', N, N, ZERO, ONE, VT, LDVT )
            CALL SLASDQ( 'U', 0, N, N, N, 0, D, E, VT, LDVT, U, LDU, U, LDU, WORK( WSTART ), INFO )
         ELSE IF( ICOMPQ.EQ.1 ) THEN
            IU = 1
            IVT = IU + N
            CALL SLASET( 'A', N, N, ZERO, ONE, Q( IU+( QSTART-1 )*N ), N )             CALL SLASET( 'A', N, N, ZERO, ONE, Q( IVT+( QSTART-1 )*N ), N )             CALL SLASDQ( 'U', 0, N, N, N, 0, D, E, Q( IVT+( QSTART-1 )*N ), N, Q( IU+( QSTART-1 )*N ), N, Q( IU+( QSTART-1 )*N ), N, WORK( WSTART ), INFO )
         END IF
         GO TO 40
      END IF

      IF( ICOMPQ.EQ.2 ) THEN
         CALL SLASET( 'A', N, N, ZERO, ONE, U, LDU )
         CALL SLASET( 'A', N, N, ZERO, ONE, VT, LDVT )
      END IF

      // Scale.

      ORGNRM = SLANST( 'M', N, D, E )
      IF( ORGNRM.EQ.ZERO ) RETURN
      CALL SLASCL( 'G', 0, 0, ORGNRM, ONE, N, 1, D, N, IERR )
      CALL SLASCL( 'G', 0, 0, ORGNRM, ONE, NM1, 1, E, NM1, IERR )

      EPS = SLAMCH( 'Epsilon' )

      MLVL = INT( LOG( REAL( N ) / REAL( SMLSIZ+1 ) ) / LOG( TWO ) ) + 1
      SMLSZP = SMLSIZ + 1

      IF( ICOMPQ.EQ.1 ) THEN
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
      END IF

      DO 20 I = 1, N
         IF( ABS( D( I ) ).LT.EPS ) THEN
            D( I ) = SIGN( EPS, D( I ) )
         END IF
   20 CONTINUE

      START = 1
      SQRE = 0

      DO 30 I = 1, NM1
         IF( ( ABS( E( I ) ).LT.EPS ) .OR. ( I.EQ.NM1 ) ) THEN

         // Subproblem found. First determine its size and then
         // apply divide and conquer on it.

            IF( I.LT.NM1 ) THEN

         // A subproblem with E(I) small for I < NM1.

               NSIZE = I - START + 1
            ELSE IF( ABS( E( I ) ).GE.EPS ) THEN

         // A subproblem with E(NM1) not too small but I = NM1.

               NSIZE = N - START + 1
            ELSE

         // A subproblem with E(NM1) small. This implies an
         // 1-by-1 subproblem at D(N). Solve this 1-by-1 problem
         // first.

               NSIZE = I - START + 1
               IF( ICOMPQ.EQ.2 ) THEN
                  U( N, N ) = SIGN( ONE, D( N ) )
                  VT( N, N ) = ONE
               ELSE IF( ICOMPQ.EQ.1 ) THEN
                  Q( N+( QSTART-1 )*N ) = SIGN( ONE, D( N ) )
                  Q( N+( SMLSIZ+QSTART-1 )*N ) = ONE
               END IF
               D( N ) = ABS( D( N ) )
            END IF
            IF( ICOMPQ.EQ.2 ) THEN
               CALL SLASD0( NSIZE, SQRE, D( START ), E( START ), U( START, START ), LDU, VT( START, START ), LDVT, SMLSIZ, IWORK, WORK( WSTART ), INFO )
            ELSE
               CALL SLASDA( ICOMPQ, SMLSIZ, NSIZE, SQRE, D( START ), E( START ), Q( START+( IU+QSTART-2 )*N ), N, Q( START+( IVT+QSTART-2 )*N ), IQ( START+K*N ), Q( START+( DIFL+QSTART-2 )* N ), Q( START+( DIFR+QSTART-2 )*N ), Q( START+( Z+QSTART-2 )*N ), Q( START+( POLES+QSTART-2 )*N ), IQ( START+GIVPTR*N ), IQ( START+GIVCOL*N ), N, IQ( START+PERM*N ), Q( START+( GIVNUM+QSTART-2 )*N ), Q( START+( IC+QSTART-2 )*N ), Q( START+( IS+QSTART-2 )*N ), WORK( WSTART ), IWORK, INFO )
            END IF
            IF( INFO.NE.0 ) THEN
               RETURN
            END IF
            START = I + 1
         END IF
   30 CONTINUE

      // Unscale

      CALL SLASCL( 'G', 0, 0, ONE, ORGNRM, N, 1, D, N, IERR )
   40 CONTINUE

      // Use Selection Sort to minimize swaps of singular vectors

      DO 60 II = 2, N
         I = II - 1
         KK = I
         P = D( I )
         DO 50 J = II, N
            IF( D( J ).GT.P ) THEN
               KK = J
               P = D( J )
            END IF
   50    CONTINUE
         IF( KK.NE.I ) THEN
            D( KK ) = D( I )
            D( I ) = P
            IF( ICOMPQ.EQ.1 ) THEN
               IQ( I ) = KK
            ELSE IF( ICOMPQ.EQ.2 ) THEN
               CALL SSWAP( N, U( 1, I ), 1, U( 1, KK ), 1 )
               CALL SSWAP( N, VT( I, 1 ), LDVT, VT( KK, 1 ), LDVT )
            END IF
         ELSE IF( ICOMPQ.EQ.1 ) THEN
            IQ( I ) = I
         END IF
   60 CONTINUE

      // If ICOMPQ = 1, use IQ(N,1) as the indicator for UPLO

      IF( ICOMPQ.EQ.1 ) THEN
         IF( IUPLO.EQ.1 ) THEN
            IQ( N ) = 1
         ELSE
            IQ( N ) = 0
         END IF
      END IF

      // If B is lower bidiagonal, update U by those Givens rotations
      // which rotated B to be upper bidiagonal

      IF( ( IUPLO.EQ.2 ) .AND. ( ICOMPQ.EQ.2 ) ) CALL SLASR( 'L', 'V', 'B', N, N, WORK( 1 ), WORK( N ), U, LDU )

      RETURN

      // End of SBDSDC

      }
