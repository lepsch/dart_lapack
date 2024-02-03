      SUBROUTINE SCHKSB( NSIZES, NN, NWDTHS, KK, NTYPES, DOTYPE, ISEED, THRESH, NOUNIT, A, LDA, SD, SE, U, LDU, WORK, LWORK, RESULT, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDU, LWORK, NOUNIT, NSIZES, NTYPES, NWDTHS;
      REAL               THRESH
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                ISEED( 4 ), KK( * ), NN( * );
      REAL               A( LDA, * ), RESULT( * ), SD( * ), SE( * ), U( LDU, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TWO, TEN
      const              ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0, TEN = 10.0E0 ;
      REAL               HALF
      const              HALF = ONE / TWO ;
      int                MAXTYP;
      const              MAXTYP = 15 ;
      // ..
      // .. Local Scalars ..
      bool               BADNN, BADNNB;
      int                I, IINFO, IMODE, ITYPE, J, JC, JCOL, JR, JSIZE, JTYPE, JWIDTH, K, KMAX, MTYPES, N, NERRS, NMATS, NMAX, NTEST, NTESTT;
      REAL               ANINV, ANORM, COND, OVFL, RTOVFL, RTUNFL, TEMP1, ULP, ULPINV, UNFL
      // ..
      // .. Local Arrays ..
      int                IDUMMA( 1 ), IOLDSD( 4 ), KMAGN( MAXTYP ), KMODE( MAXTYP ), KTYPE( MAXTYP );
      // ..
      // .. External Functions ..
      REAL               SLAMCH
      // EXTERNAL SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLACPY, SLASUM, SLATMR, SLATMS, SLASET, SSBT21, SSBTRD, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, REAL, SQRT
      // ..
      // .. Data statements ..
      DATA               KTYPE / 1, 2, 5*4, 5*5, 3*8 /
      DATA               KMAGN / 2*1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3 /       DATA               KMODE / 2*0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0 /
      // ..
      // .. Executable Statements ..

      // Check for errors

      NTESTT = 0
      INFO = 0

      // Important constants

      BADNN = .FALSE.
      NMAX = 1
      DO 10 J = 1, NSIZES
         NMAX = MAX( NMAX, NN( J ) )
         IF( NN( J ).LT.0 ) BADNN = .TRUE.
   10 CONTINUE

      BADNNB = .FALSE.
      KMAX = 0
      DO 20 J = 1, NSIZES
         KMAX = MAX( KMAX, KK( J ) )
         IF( KK( J ).LT.0 ) BADNNB = .TRUE.
   20 CONTINUE
      KMAX = MIN( NMAX-1, KMAX )

      // Check for errors

      if ( NSIZES.LT.0 ) {
         INFO = -1
      } else if ( BADNN ) {
         INFO = -2
      } else if ( NWDTHS.LT.0 ) {
         INFO = -3
      } else if ( BADNNB ) {
         INFO = -4
      } else if ( NTYPES.LT.0 ) {
         INFO = -5
      } else if ( LDA.LT.KMAX+1 ) {
         INFO = -11
      } else if ( LDU.LT.NMAX ) {
         INFO = -15
      } else if ( ( MAX( LDA, NMAX )+1 )*NMAX.GT.LWORK ) {
         INFO = -17
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'SCHKSB', -INFO )
         RETURN
      }

      // Quick return if possible

      IF( NSIZES.EQ.0 .OR. NTYPES.EQ.0 .OR. NWDTHS.EQ.0 ) RETURN

      // More Important constants

      UNFL = SLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )
      ULPINV = ONE / ULP
      RTUNFL = SQRT( UNFL )
      RTOVFL = SQRT( OVFL )

      // Loop over sizes, types

      NERRS = 0
      NMATS = 0

      DO 190 JSIZE = 1, NSIZES
         N = NN( JSIZE )
         ANINV = ONE / REAL( MAX( 1, N ) )

         DO 180 JWIDTH = 1, NWDTHS
            K = KK( JWIDTH )
            IF( K.GT.N ) GO TO 180
            K = MAX( 0, MIN( N-1, K ) )

            if ( NSIZES.NE.1 ) {
               MTYPES = MIN( MAXTYP, NTYPES )
            } else {
               MTYPES = MIN( MAXTYP+1, NTYPES )
            }

            DO 170 JTYPE = 1, MTYPES
               IF( .NOT.DOTYPE( JTYPE ) ) GO TO 170
               NMATS = NMATS + 1
               NTEST = 0

               DO 30 J = 1, 4
                  IOLDSD( J ) = ISEED( J )
   30          CONTINUE

               // Compute "A".
               // Store as "Upper"; later, we will copy to other format.

               // Control parameters:

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

               IF( MTYPES.GT.MAXTYP ) GO TO 100

               ITYPE = KTYPE( JTYPE )
               IMODE = KMODE( JTYPE )

               // Compute norm

               GO TO ( 40, 50, 60 )KMAGN( JTYPE )

   40          CONTINUE
               ANORM = ONE
               GO TO 70

   50          CONTINUE
               ANORM = ( RTOVFL*ULP )*ANINV
               GO TO 70

   60          CONTINUE
               ANORM = RTUNFL*N*ULPINV
               GO TO 70

   70          CONTINUE

               CALL SLASET( 'Full', LDA, N, ZERO, ZERO, A, LDA )
               IINFO = 0
               if ( JTYPE.LE.15 ) {
                  COND = ULPINV
               } else {
                  COND = ULPINV*ANINV / TEN
               }

               // Special Matrices -- Identity & Jordan block

                  // Zero

               if ( ITYPE.EQ.1 ) {
                  IINFO = 0

               } else if ( ITYPE.EQ.2 ) {

                  // Identity

                  DO 80 JCOL = 1, N
                     A( K+1, JCOL ) = ANORM
   80             CONTINUE

               } else if ( ITYPE.EQ.4 ) {

                  // Diagonal Matrix, [Eigen]values Specified

                  CALL SLATMS( N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, 0, 0, 'Q', A( K+1, 1 ), LDA, WORK( N+1 ), IINFO )

               } else if ( ITYPE.EQ.5 ) {

                  // Symmetric, eigenvalues specified

                  CALL SLATMS( N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, K, K, 'Q', A, LDA, WORK( N+1 ), IINFO )

               } else if ( ITYPE.EQ.7 ) {

                  // Diagonal, random eigenvalues

                  CALL SLATMR( N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, 0, 0, ZERO, ANORM, 'Q', A( K+1, 1 ), LDA, IDUMMA, IINFO )

               } else if ( ITYPE.EQ.8 ) {

                  // Symmetric, random eigenvalues

                  CALL SLATMR( N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, K, K, ZERO, ANORM, 'Q', A, LDA, IDUMMA, IINFO )

               } else if ( ITYPE.EQ.9 ) {

                  // Positive definite, eigenvalues specified.

                  CALL SLATMS( N, N, 'S', ISEED, 'P', WORK, IMODE, COND, ANORM, K, K, 'Q', A, LDA, WORK( N+1 ), IINFO )

               } else if ( ITYPE.EQ.10 ) {

                  // Positive definite tridiagonal, eigenvalues specified.

                  IF( N.GT.1 ) K = MAX( 1, K )                   CALL SLATMS( N, N, 'S', ISEED, 'P', WORK, IMODE, COND, ANORM, 1, 1, 'Q', A( K, 1 ), LDA, WORK( N+1 ), IINFO )
                  DO 90 I = 2, N
                     TEMP1 = ABS( A( K, I ) ) / SQRT( ABS( A( K+1, I-1 )*A( K+1, I ) ) )
                     if ( TEMP1.GT.HALF ) {
                        A( K, I ) = HALF*SQRT( ABS( A( K+1, I-1 )*A( K+1, I ) ) )
                     }
   90             CONTINUE

               } else {

                  IINFO = 1
               }

               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'Generator', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               }

  100          CONTINUE

               // Call SSBTRD to compute S and U from upper triangle.

               CALL SLACPY( ' ', K+1, N, A, LDA, WORK, LDA )

               NTEST = 1
               CALL SSBTRD( 'V', 'U', N, K, WORK, LDA, SD, SE, U, LDU, WORK( LDA*N+1 ), IINFO )

               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSBTRD(U)', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( 1 ) = ULPINV
                     GO TO 150
                  }
               }

               // Do tests 1 and 2

               CALL SSBT21( 'Upper', N, K, 1, A, LDA, SD, SE, U, LDU, WORK, RESULT( 1 ) )

               // Convert A from Upper-Triangle-Only storage to
               // Lower-Triangle-Only storage.

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

               // Call SSBTRD to compute S and U from lower triangle

               CALL SLACPY( ' ', K+1, N, A, LDA, WORK, LDA )

               NTEST = 3
               CALL SSBTRD( 'V', 'L', N, K, WORK, LDA, SD, SE, U, LDU, WORK( LDA*N+1 ), IINFO )

               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSBTRD(L)', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( 3 ) = ULPINV
                     GO TO 150
                  }
               }
               NTEST = 4

               // Do tests 3 and 4

               CALL SSBT21( 'Lower', N, K, 1, A, LDA, SD, SE, U, LDU, WORK, RESULT( 3 ) )

               // End of Loop -- Check for RESULT(j) > THRESH

  150          CONTINUE
               NTESTT = NTESTT + NTEST

               // Print out tests which fail.

               DO 160 JR = 1, NTEST
                  if ( RESULT( JR ).GE.THRESH ) {

                     // If this is the first test to fail,
                     // print a header to the data file.

                     if ( NERRS.EQ.0 ) {
                        WRITE( NOUNIT, FMT = 9998 )'SSB'
                        WRITE( NOUNIT, FMT = 9997 )
                        WRITE( NOUNIT, FMT = 9996 )
                        WRITE( NOUNIT, FMT = 9995 )'Symmetric'
                        WRITE( NOUNIT, FMT = 9994 )'orthogonal', '''', 'transpose', ( '''', J = 1, 4 )
                     }
                     NERRS = NERRS + 1
                     WRITE( NOUNIT, FMT = 9993 )N, K, IOLDSD, JTYPE, JR, RESULT( JR )
                  }
  160          CONTINUE

  170       CONTINUE
  180    CONTINUE
  190 CONTINUE

      // Summary

      CALL SLASUM( 'SSB', NOUNIT, NERRS, NTESTT )
      RETURN

 9999 FORMAT( ' SCHKSB: ', A, ' returned INFO=', I6, '.', / 9X, 'N=',
     $      I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )

 9998 FORMAT( / 1X, A3,
     $      ' -- Real Symmetric Banded Tridiagonal Reduction Routines' )
 9997 FORMAT( ' Matrix types (see SCHKSB for details): ' )

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

 9994 FORMAT( / ' Tests performed:   (S is Tridiag,  U is ', A, ',',
     $      / 20X, A, ' means ', A, '.', / ' UPLO=''U'':',
     $      / '  1= | A - U S U', A1, ' | / ( |A| n ulp )     ',
     $      '  2= | I - U U', A1, ' | / ( n ulp )', / ' UPLO=''L'':',
     $      / '  3= | A - U S U', A1, ' | / ( |A| n ulp )     ',
     $      '  4= | I - U U', A1, ' | / ( n ulp )' )
 9993 FORMAT( ' N=', I5, ', K=', I4, ', seed=', 4( I4, ',' ), ' type ',
     $      I2, ', test(', I2, ')=', G10.3 )

      // End of SCHKSB

      }
