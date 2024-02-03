      SUBROUTINE DCHKSB( NSIZES, NN, NWDTHS, KK, NTYPES, DOTYPE, ISEED, THRESH, NOUNIT, A, LDA, SD, SE, U, LDU, WORK, LWORK, RESULT, INFO )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                INFO, LDA, LDU, LWORK, NOUNIT, NSIZES, NTYPES, NWDTHS;
      double             THRESH;
*     ..
*     .. Array Arguments ..
      bool               DOTYPE( * );
      int                ISEED( 4 ), KK( * ), NN( * );
      double             A( LDA, * ), RESULT( * ), SD( * ), SE( * ), U( LDU, * ), WORK( * );
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double             ZERO, ONE, TWO, TEN;
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, TEN = 10.0D0 )
      double             HALF;
      PARAMETER          ( HALF = ONE / TWO )
      int                MAXTYP;
      PARAMETER          ( MAXTYP = 15 )
*     ..
*     .. Local Scalars ..
      bool               BADNN, BADNNB;
      int                I, IINFO, IMODE, ITYPE, J, JC, JCOL, JR, JSIZE, JTYPE, JWIDTH, K, KMAX, MTYPES, N, NERRS, NMATS, NMAX, NTEST, NTESTT;
      double             ANINV, ANORM, COND, OVFL, RTOVFL, RTUNFL, TEMP1, ULP, ULPINV, UNFL;
*     ..
*     .. Local Arrays ..
      int                IDUMMA( 1 ), IOLDSD( 4 ), KMAGN( MAXTYP ), KMODE( MAXTYP ), KTYPE( MAXTYP );
*     ..
*     .. External Functions ..
      double             DLAMCH;
      EXTERNAL           DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLACPY, DLASET, DLASUM, DLATMR, DLATMS, DSBT21, DSBTRD, XERBLA
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, MIN, SQRT
*     ..
*     .. Data statements ..
      DATA               KTYPE / 1, 2, 5*4, 5*5, 3*8 /
      DATA               KMAGN / 2*1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3 /       DATA               KMODE / 2*0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0 /
*     ..
*     .. Executable Statements ..
*
*     Check for errors
*
      NTESTT = 0
      INFO = 0
*
*     Important constants
*
      BADNN = .FALSE.
      NMAX = 1
      DO 10 J = 1, NSIZES
         NMAX = MAX( NMAX, NN( J ) )
         IF( NN( J ).LT.0 ) BADNN = .TRUE.
   10 CONTINUE
*
      BADNNB = .FALSE.
      KMAX = 0
      DO 20 J = 1, NSIZES
         KMAX = MAX( KMAX, KK( J ) )
         IF( KK( J ).LT.0 ) BADNNB = .TRUE.
   20 CONTINUE
      KMAX = MIN( NMAX-1, KMAX )
*
*     Check for errors
*
      IF( NSIZES.LT.0 ) THEN
         INFO = -1
      ELSE IF( BADNN ) THEN
         INFO = -2
      ELSE IF( NWDTHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( BADNNB ) THEN
         INFO = -4
      ELSE IF( NTYPES.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.KMAX+1 ) THEN
         INFO = -11
      ELSE IF( LDU.LT.NMAX ) THEN
         INFO = -15
      ELSE IF( ( MAX( LDA, NMAX )+1 )*NMAX.GT.LWORK ) THEN
         INFO = -17
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DCHKSB', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( NSIZES.EQ.0 .OR. NTYPES.EQ.0 .OR. NWDTHS.EQ.0 ) RETURN
*
*     More Important constants
*
      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )
      ULPINV = ONE / ULP
      RTUNFL = SQRT( UNFL )
      RTOVFL = SQRT( OVFL )
*
*     Loop over sizes, types
*
      NERRS = 0
      NMATS = 0
*
      DO 190 JSIZE = 1, NSIZES
         N = NN( JSIZE )
         ANINV = ONE / DBLE( MAX( 1, N ) )
*
         DO 180 JWIDTH = 1, NWDTHS
            K = KK( JWIDTH )
            IF( K.GT.N ) GO TO 180
            K = MAX( 0, MIN( N-1, K ) )
*
            IF( NSIZES.NE.1 ) THEN
               MTYPES = MIN( MAXTYP, NTYPES )
            ELSE
               MTYPES = MIN( MAXTYP+1, NTYPES )
            END IF
*
            DO 170 JTYPE = 1, MTYPES
               IF( .NOT.DOTYPE( JTYPE ) ) GO TO 170
               NMATS = NMATS + 1
               NTEST = 0
*
               DO 30 J = 1, 4
                  IOLDSD( J ) = ISEED( J )
   30          CONTINUE
*
*              Compute "A".
*              Store as "Upper"; later, we will copy to other format.
*
*              Control parameters:
*
*                  KMAGN  KMODE        KTYPE
*              =1  O(1)   clustered 1  zero
*              =2  large  clustered 2  identity
*              =3  small  exponential  (none)
*              =4         arithmetic   diagonal, (w/ eigenvalues)
*              =5         random log   symmetric, w/ eigenvalues
*              =6         random       (none)
*              =7                      random diagonal
*              =8                      random symmetric
*              =9                      positive definite
*              =10                     diagonally dominant tridiagonal
*
               IF( MTYPES.GT.MAXTYP ) GO TO 100
*
               ITYPE = KTYPE( JTYPE )
               IMODE = KMODE( JTYPE )
*
*              Compute norm
*
               GO TO ( 40, 50, 60 )KMAGN( JTYPE )
*
   40          CONTINUE
               ANORM = ONE
               GO TO 70
*
   50          CONTINUE
               ANORM = ( RTOVFL*ULP )*ANINV
               GO TO 70
*
   60          CONTINUE
               ANORM = RTUNFL*N*ULPINV
               GO TO 70
*
   70          CONTINUE
*
               CALL DLASET( 'Full', LDA, N, ZERO, ZERO, A, LDA )
               IINFO = 0
               IF( JTYPE.LE.15 ) THEN
                  COND = ULPINV
               ELSE
                  COND = ULPINV*ANINV / TEN
               END IF
*
*              Special Matrices -- Identity & Jordan block
*
*                 Zero
*
               IF( ITYPE.EQ.1 ) THEN
                  IINFO = 0
*
               ELSE IF( ITYPE.EQ.2 ) THEN
*
*                 Identity
*
                  DO 80 JCOL = 1, N
                     A( K+1, JCOL ) = ANORM
   80             CONTINUE
*
               ELSE IF( ITYPE.EQ.4 ) THEN
*
*                 Diagonal Matrix, [Eigen]values Specified
*
                  CALL DLATMS( N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, 0, 0, 'Q', A( K+1, 1 ), LDA, WORK( N+1 ), IINFO )
*
               ELSE IF( ITYPE.EQ.5 ) THEN
*
*                 Symmetric, eigenvalues specified
*
                  CALL DLATMS( N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, K, K, 'Q', A, LDA, WORK( N+1 ), IINFO )
*
               ELSE IF( ITYPE.EQ.7 ) THEN
*
*                 Diagonal, random eigenvalues
*
                  CALL DLATMR( N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, 0, 0, ZERO, ANORM, 'Q', A( K+1, 1 ), LDA, IDUMMA, IINFO )
*
               ELSE IF( ITYPE.EQ.8 ) THEN
*
*                 Symmetric, random eigenvalues
*
                  CALL DLATMR( N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, K, K, ZERO, ANORM, 'Q', A, LDA, IDUMMA, IINFO )
*
               ELSE IF( ITYPE.EQ.9 ) THEN
*
*                 Positive definite, eigenvalues specified.
*
                  CALL DLATMS( N, N, 'S', ISEED, 'P', WORK, IMODE, COND, ANORM, K, K, 'Q', A, LDA, WORK( N+1 ), IINFO )
*
               ELSE IF( ITYPE.EQ.10 ) THEN
*
*                 Positive definite tridiagonal, eigenvalues specified.
*
                  IF( N.GT.1 ) K = MAX( 1, K )                   CALL DLATMS( N, N, 'S', ISEED, 'P', WORK, IMODE, COND, ANORM, 1, 1, 'Q', A( K, 1 ), LDA, WORK( N+1 ), IINFO )
                  DO 90 I = 2, N
                     TEMP1 = ABS( A( K, I ) ) / SQRT( ABS( A( K+1, I-1 )*A( K+1, I ) ) )
                     IF( TEMP1.GT.HALF ) THEN
                        A( K, I ) = HALF*SQRT( ABS( A( K+1, I-1 )*A( K+1, I ) ) )
                     END IF
   90             CONTINUE
*
               ELSE
*
                  IINFO = 1
               END IF
*
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'Generator', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               END IF
*
  100          CONTINUE
*
*              Call DSBTRD to compute S and U from upper triangle.
*
               CALL DLACPY( ' ', K+1, N, A, LDA, WORK, LDA )
*
               NTEST = 1
               CALL DSBTRD( 'V', 'U', N, K, WORK, LDA, SD, SE, U, LDU, WORK( LDA*N+1 ), IINFO )
*
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSBTRD(U)', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( 1 ) = ULPINV
                     GO TO 150
                  END IF
               END IF
*
*              Do tests 1 and 2
*
               CALL DSBT21( 'Upper', N, K, 1, A, LDA, SD, SE, U, LDU, WORK, RESULT( 1 ) )
*
*              Convert A from Upper-Triangle-Only storage to
*              Lower-Triangle-Only storage.
*
               DO 120 JC = 1, N
                  DO 110 JR = 0, MIN( K, N-JC )
                     A( JR+1, JC ) = A( K+1-JR, JC+JR )
  110             CONTINUE
  120          CONTINUE
               DO 140 JC = N + 1 - K, N
                  DO 130 JR = MIN( K, N-JC ) + 1, K
                     A( JR+1, JC ) = ZERO
  130             CONTINUE
  140          CONTINUE
*
*              Call DSBTRD to compute S and U from lower triangle
*
               CALL DLACPY( ' ', K+1, N, A, LDA, WORK, LDA )
*
               NTEST = 3
               CALL DSBTRD( 'V', 'L', N, K, WORK, LDA, SD, SE, U, LDU, WORK( LDA*N+1 ), IINFO )
*
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'DSBTRD(L)', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( 3 ) = ULPINV
                     GO TO 150
                  END IF
               END IF
               NTEST = 4
*
*              Do tests 3 and 4
*
               CALL DSBT21( 'Lower', N, K, 1, A, LDA, SD, SE, U, LDU, WORK, RESULT( 3 ) )
*
*              End of Loop -- Check for RESULT(j) > THRESH
*
  150          CONTINUE
               NTESTT = NTESTT + NTEST
*
*              Print out tests which fail.
*
               DO 160 JR = 1, NTEST
                  IF( RESULT( JR ).GE.THRESH ) THEN
*
*                    If this is the first test to fail,
*                    print a header to the data file.
*
                     IF( NERRS.EQ.0 ) THEN
                        WRITE( NOUNIT, FMT = 9998 )'DSB'
                        WRITE( NOUNIT, FMT = 9997 )
                        WRITE( NOUNIT, FMT = 9996 )
                        WRITE( NOUNIT, FMT = 9995 )'Symmetric'
                        WRITE( NOUNIT, FMT = 9994 )'orthogonal', '''', 'transpose', ( '''', J = 1, 4 )
                     END IF
                     NERRS = NERRS + 1
                     WRITE( NOUNIT, FMT = 9993 )N, K, IOLDSD, JTYPE, JR, RESULT( JR )
                  END IF
  160          CONTINUE
*
  170       CONTINUE
  180    CONTINUE
  190 CONTINUE
*
*     Summary
*
      CALL DLASUM( 'DSB', NOUNIT, NERRS, NTESTT )
      RETURN
*
 9999 FORMAT( ' DCHKSB: ', A, ' returned INFO=', I6, '.', / 9X, 'N=',
     $      I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )
*
 9998 FORMAT( / 1X, A3,
     $      ' -- Real Symmetric Banded Tridiagonal Reduction Routines' )
 9997 FORMAT( ' Matrix types (see DCHKSB for details): ' )
*
 9996 FORMAT( / ' Special Matrices:',
     $      / '  1=Zero matrix.                        ',
     $      '  5=Diagonal: clustered entries.',
     $      / '  2=Identity matrix.                    ',
     $      '  6=Diagonal: large, evenly spaced.',
     $      / '  3=Diagonal: evenly spaced entries.    ',
     $      '  7=Diagonal: small, evenly spaced.',
     $      / '  4=Diagonal: geometr. spaced entries.' )
 9995 FORMAT( ' Dense ', A, ' Banded Matrices:',
     $      / '  8=Evenly spaced eigenvals.            ',
     $      ' 12=Small, evenly spaced eigenvals.',
     $      / '  9=Geometrically spaced eigenvals.     ',
     $      ' 13=Matrix with random O(1) entries.',
     $      / ' 10=Clustered eigenvalues.              ',
     $      ' 14=Matrix with large random entries.',
     $      / ' 11=Large, evenly spaced eigenvals.     ',
     $      ' 15=Matrix with small random entries.' )
*
 9994 FORMAT( / ' Tests performed:   (S is Tridiag,  U is ', A, ',',
     $      / 20X, A, ' means ', A, '.', / ' UPLO=''U'':',
     $      / '  1= | A - U S U', A1, ' | / ( |A| n ulp )     ',
     $      '  2= | I - U U', A1, ' | / ( n ulp )', / ' UPLO=''L'':',
     $      / '  3= | A - U S U', A1, ' | / ( |A| n ulp )     ',
     $      '  4= | I - U U', A1, ' | / ( n ulp )' )
 9993 FORMAT( ' N=', I5, ', K=', I4, ', seed=', 4( I4, ',' ), ' type ',
     $      I2, ', test(', I2, ')=', G10.3 )
*
*     End of DCHKSB
*
      END
