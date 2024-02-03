      SUBROUTINE SCHKSB2STG( NSIZES, NN, NWDTHS, KK, NTYPES, DOTYPE, ISEED, THRESH, NOUNIT, A, LDA, SD, SE, D1, D2, D3, U, LDU, WORK, LWORK, RESULT, INFO )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      // .. Scalar Arguments ..
      int                INFO, LDA, LDU, LWORK, NOUNIT, NSIZES, NTYPES, NWDTHS;
      REAL               THRESH
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                ISEED( 4 ), KK( * ), NN( * );
      REAL               A( LDA, * ), RESULT( * ), SD( * ), SE( * ), D1( * ), D2( * ), D3( * ), U( LDU, * ), WORK( * )
      // ..
*
*  =====================================================================
*
      // .. Parameters ..
      REAL               ZERO, ONE, TWO, TEN
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0, TEN = 10.0E0 )
      REAL               HALF
      PARAMETER          ( HALF = ONE / TWO )
      int                MAXTYP;
      PARAMETER          ( MAXTYP = 15 )
      // ..
      // .. Local Scalars ..
      bool               BADNN, BADNNB;
      int                I, IINFO, IMODE, ITYPE, J, JC, JCOL, JR, JSIZE, JTYPE, JWIDTH, K, KMAX, LH, LW, MTYPES, N, NERRS, NMATS, NMAX, NTEST, NTESTT;
      REAL               ANINV, ANORM, COND, OVFL, RTOVFL, RTUNFL, TEMP1, TEMP2, TEMP3, TEMP4, ULP, ULPINV, UNFL
      // ..
      // .. Local Arrays ..
      int                IDUMMA( 1 ), IOLDSD( 4 ), KMAGN( MAXTYP ), KMODE( MAXTYP ), KTYPE( MAXTYP );
      // ..
      // .. External Functions ..
      REAL               SLAMCH
      // EXTERNAL SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLACPY, SLASET, SLASUM, SLATMR, SLATMS, SSBT21, SSBTRD, XERBLA, SSYTRD_SB2ST, SSTEQR
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, REAL, MAX, MIN, SQRT
      // ..
      // .. Data statements ..
      DATA               KTYPE / 1, 2, 5*4, 5*5, 3*8 /
      DATA               KMAGN / 2*1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3 /       DATA               KMODE / 2*0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0 /
      // ..
      // .. Executable Statements ..
*
      // Check for errors
*
      NTESTT = 0
      INFO = 0
*
      // Important constants
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
      // Check for errors
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
         CALL XERBLA( 'SCHKSB2STG', -INFO )
         RETURN
      END IF
*
      // Quick return if possible
*
      IF( NSIZES.EQ.0 .OR. NTYPES.EQ.0 .OR. NWDTHS.EQ.0 ) RETURN
*
      // More Important constants
*
      UNFL = SLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )
      ULPINV = ONE / ULP
      RTUNFL = SQRT( UNFL )
      RTOVFL = SQRT( OVFL )
*
      // Loop over sizes, types
*
      NERRS = 0
      NMATS = 0
*
      DO 190 JSIZE = 1, NSIZES
         N = NN( JSIZE )
         ANINV = ONE / REAL( MAX( 1, N ) )
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
               // Compute "A".
               // Store as "Upper"; later, we will copy to other format.
*
               // Control parameters:
*
                   // KMAGN  KMODE        KTYPE
               // =1  O(1)   clustered 1  zero
               // =2  large  clustered 2  identity
               // =3  small  exponential  (none)
               // =4         arithmetic   diagonal, (w/ eigenvalues)
               // =5         random log   symmetric, w/ eigenvalues
               // =6         random       (none)
               // =7                      random diagonal
               // =8                      random symmetric
               // =9                      positive definite
               // =10                     diagonally dominant tridiagonal
*
               IF( MTYPES.GT.MAXTYP ) GO TO 100
*
               ITYPE = KTYPE( JTYPE )
               IMODE = KMODE( JTYPE )
*
               // Compute norm
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
               CALL SLASET( 'Full', LDA, N, ZERO, ZERO, A, LDA )
               IINFO = 0
               IF( JTYPE.LE.15 ) THEN
                  COND = ULPINV
               ELSE
                  COND = ULPINV*ANINV / TEN
               END IF
*
               // Special Matrices -- Identity & Jordan block
*
                  // Zero
*
               IF( ITYPE.EQ.1 ) THEN
                  IINFO = 0
*
               ELSE IF( ITYPE.EQ.2 ) THEN
*
                  // Identity
*
                  DO 80 JCOL = 1, N
                     A( K+1, JCOL ) = ANORM
   80             CONTINUE
*
               ELSE IF( ITYPE.EQ.4 ) THEN
*
                  // Diagonal Matrix, [Eigen]values Specified
*
                  CALL SLATMS( N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, 0, 0, 'Q', A( K+1, 1 ), LDA, WORK( N+1 ), IINFO )
*
               ELSE IF( ITYPE.EQ.5 ) THEN
*
                  // Symmetric, eigenvalues specified
*
                  CALL SLATMS( N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, K, K, 'Q', A, LDA, WORK( N+1 ), IINFO )
*
               ELSE IF( ITYPE.EQ.7 ) THEN
*
                  // Diagonal, random eigenvalues
*
                  CALL SLATMR( N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, 0, 0, ZERO, ANORM, 'Q', A( K+1, 1 ), LDA, IDUMMA, IINFO )
*
               ELSE IF( ITYPE.EQ.8 ) THEN
*
                  // Symmetric, random eigenvalues
*
                  CALL SLATMR( N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, K, K, ZERO, ANORM, 'Q', A, LDA, IDUMMA, IINFO )
*
               ELSE IF( ITYPE.EQ.9 ) THEN
*
                  // Positive definite, eigenvalues specified.
*
                  CALL SLATMS( N, N, 'S', ISEED, 'P', WORK, IMODE, COND, ANORM, K, K, 'Q', A, LDA, WORK( N+1 ), IINFO )
*
               ELSE IF( ITYPE.EQ.10 ) THEN
*
                  // Positive definite tridiagonal, eigenvalues specified.
*
                  IF( N.GT.1 ) K = MAX( 1, K )                   CALL SLATMS( N, N, 'S', ISEED, 'P', WORK, IMODE, COND, ANORM, 1, 1, 'Q', A( K, 1 ), LDA, WORK( N+1 ), IINFO )
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
               // Call SSBTRD to compute S and U from upper triangle.
*
               CALL SLACPY( ' ', K+1, N, A, LDA, WORK, LDA )
*
               NTEST = 1
               CALL SSBTRD( 'V', 'U', N, K, WORK, LDA, SD, SE, U, LDU, WORK( LDA*N+1 ), IINFO )
*
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSBTRD(U)', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( 1 ) = ULPINV
                     GO TO 150
                  END IF
               END IF
*
               // Do tests 1 and 2
*
               CALL SSBT21( 'Upper', N, K, 1, A, LDA, SD, SE, U, LDU, WORK, RESULT( 1 ) )
*
               // Before converting A into lower for SSBTRD, run SSYTRD_SB2ST
               // otherwise matrix A will be converted to lower and then need
              t // o be converted back to upper in order to run the upper case
               // ofSSYTRD_SB2ST
*
               // Compute D1 the eigenvalues resulting from the tridiagonal
               // form using the SSBTRD and used as reference to compare
               // with the SSYTRD_SB2ST routine
*
               // Compute D1 from the SSBTRD and used as reference for the
               // SSYTRD_SB2ST
*
               CALL SCOPY( N, SD, 1, D1, 1 )
               IF( N.GT.0 ) CALL SCOPY( N-1, SE, 1, WORK, 1 )
*
               CALL SSTEQR( 'N', N, D1, WORK, WORK( N+1 ), LDU, WORK( N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSTEQR(N)', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( 5 ) = ULPINV
                     GO TO 150
                  END IF
               END IF
*
               // SSYTRD_SB2ST Upper case is used to compute D2.
               // Note to set SD and SE to zero to be sure not reusing
              t // he one from above. Compare it with D1 computed
               // using the SSBTRD.
*
               CALL SLASET( 'Full', N, 1, ZERO, ZERO, SD, N )
               CALL SLASET( 'Full', N, 1, ZERO, ZERO, SE, N )
               CALL SLACPY( ' ', K+1, N, A, LDA, U, LDU )
               LH = MAX(1, 4*N)
               LW = LWORK - LH
               CALL SSYTRD_SB2ST( 'N', 'N', "U", N, K, U, LDU, SD, SE, WORK, LH, WORK( LH+1 ), LW, IINFO )
*
               // Compute D2 from the SSYTRD_SB2ST Upper case
*
               CALL SCOPY( N, SD, 1, D2, 1 )
               IF( N.GT.0 ) CALL SCOPY( N-1, SE, 1, WORK, 1 )
*
               CALL SSTEQR( 'N', N, D2, WORK, WORK( N+1 ), LDU, WORK( N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSTEQR(N)', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( 5 ) = ULPINV
                     GO TO 150
                  END IF
               END IF
*
               // Convert A from Upper-Triangle-Only storage to
               // Lower-Triangle-Only storage.
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
               // Call SSBTRD to compute S and U from lower triangle
*
               CALL SLACPY( ' ', K+1, N, A, LDA, WORK, LDA )
*
               NTEST = 3
               CALL SSBTRD( 'V', 'L', N, K, WORK, LDA, SD, SE, U, LDU, WORK( LDA*N+1 ), IINFO )
*
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSBTRD(L)', IINFO, N, JTYPE, IOLDSD
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
               // Do tests 3 and 4
*
               CALL SSBT21( 'Lower', N, K, 1, A, LDA, SD, SE, U, LDU, WORK, RESULT( 3 ) )
*
               // SSYTRD_SB2ST Lower case is used to compute D3.
               // Note to set SD and SE to zero to be sure not reusing
              t // he one from above. Compare it with D1 computed
               // using the SSBTRD.
*
               CALL SLASET( 'Full', N, 1, ZERO, ZERO, SD, N )
               CALL SLASET( 'Full', N, 1, ZERO, ZERO, SE, N )
               CALL SLACPY( ' ', K+1, N, A, LDA, U, LDU )
               LH = MAX(1, 4*N)
               LW = LWORK - LH
               CALL SSYTRD_SB2ST( 'N', 'N', "L", N, K, U, LDU, SD, SE, WORK, LH, WORK( LH+1 ), LW, IINFO )
*
               // Compute D3 from the 2-stage Upper case
*
               CALL SCOPY( N, SD, 1, D3, 1 )
               IF( N.GT.0 ) CALL SCOPY( N-1, SE, 1, WORK, 1 )
*
               CALL SSTEQR( 'N', N, D3, WORK, WORK( N+1 ), LDU, WORK( N+1 ), IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )'SSTEQR(N)', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  IF( IINFO.LT.0 ) THEN
                     RETURN
                  ELSE
                     RESULT( 6 ) = ULPINV
                     GO TO 150
                  END IF
               END IF
*
*
               // Do Tests 3 and 4 which are similar to 11 and 12 but with the
               // D1 computed using the standard 1-stage reduction as reference
*
               NTEST = 6
               TEMP1 = ZERO
               TEMP2 = ZERO
               TEMP3 = ZERO
               TEMP4 = ZERO
*
               DO 151 J = 1, N
                  TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D2( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( D1( J )-D2( J ) ) )
                  TEMP3 = MAX( TEMP3, ABS( D1( J ) ), ABS( D3( J ) ) )
                  TEMP4 = MAX( TEMP4, ABS( D1( J )-D3( J ) ) )
  151          CONTINUE
*
               RESULT(5) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )
               RESULT(6) = TEMP4 / MAX( UNFL, ULP*MAX( TEMP3, TEMP4 ) )
*
               // End of Loop -- Check for RESULT(j) > THRESH
*
  150          CONTINUE
               NTESTT = NTESTT + NTEST
*
               // Print out tests which fail.
*
               DO 160 JR = 1, NTEST
                  IF( RESULT( JR ).GE.THRESH ) THEN
*
                     // If this is the first test to fail,
                     // print a header to the data file.
*
                     IF( NERRS.EQ.0 ) THEN
                        WRITE( NOUNIT, FMT = 9998 )'SSB'
                        WRITE( NOUNIT, FMT = 9997 )
                        WRITE( NOUNIT, FMT = 9996 )
                        WRITE( NOUNIT, FMT = 9995 )'Symmetric'
                        WRITE( NOUNIT, FMT = 9994 )'orthogonal', '''', 'transpose', ( '''', J = 1, 6 )
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
      // Summary
*
      CALL SLASUM( 'SSB', NOUNIT, NERRS, NTESTT )
      RETURN
*
 9999 FORMAT( ' SCHKSB2STG: ', A, ' returned INFO=', I6, '.', / 9X,
     $      'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5,
     $      ')' )
*
 9998 FORMAT( / 1X, A3,
     $      ' -- Real Symmetric Banded Tridiagonal Reduction Routines' )
 9997 FORMAT( ' Matrix types (see SCHKSB2STG for details): ' )
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
     $      '  4= | I - U U', A1, ' | / ( n ulp )' / ' Eig check:',
     $      /'  5= | D1 - D2', '', ' | / ( |D1| ulp )         ',
     $      '  6= | D1 - D3', '', ' | / ( |D1| ulp )          ' )
 9993 FORMAT( ' N=', I5, ', K=', I4, ', seed=', 4( I4, ',' ), ' type ',
     $      I2, ', test(', I2, ')=', G10.3 )
*
      // End of SCHKSB2STG
*
      END
