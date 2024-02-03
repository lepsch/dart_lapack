      SUBROUTINE DDRVVX( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH, NIUNIT, NOUNIT, A, LDA, H, WR, WI, WR1, WI1, VL, LDVL, VR, LDVR, LRE, LDLRE, RCONDV, RCNDV1, RCDVIN, RCONDE, RCNDE1, RCDEIN, SCALE, SCALE1, RESULT, WORK, NWORK, IWORK, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDLRE, LDVL, LDVR, NIUNIT, NOUNIT, NSIZES, NTYPES, NWORK;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                ISEED( 4 ), IWORK( * ), NN( * );
      double             A( LDA, * ), H( LDA, * ), LRE( LDLRE, * ), RCDEIN( * ), RCDVIN( * ), RCNDE1( * ), RCNDV1( * ), RCONDE( * ), RCONDV( * ), RESULT( 11 ), SCALE( * ), SCALE1( * ), VL( LDVL, * ), VR( LDVR, * ), WI( * ), WI1( * ), WORK( * ), WR( * ), WR1( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
      int                MAXTYP;
      const              MAXTYP = 21 ;
      // ..
      // .. Local Scalars ..
      bool               BADNN;
      String             BALANC;
      String             PATH;
      int                I, IBAL, IINFO, IMODE, ITYPE, IWK, J, JCOL, JSIZE, JTYPE, MTYPES, N, NERRS, NFAIL, NMAX, NNWORK, NTEST, NTESTF, NTESTT;
      double             ANORM, COND, CONDS, OVFL, RTULP, RTULPI, ULP, ULPINV, UNFL;
      // ..
      // .. Local Arrays ..
      String             ADUMMA( 1 ), BAL( 4 );
      int                IDUMMA( 1 ), IOLDSD( 4 ), KCONDS( MAXTYP ), KMAGN( MAXTYP ), KMODE( MAXTYP ), KTYPE( MAXTYP );
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DGET23, DLASET, DLASUM, DLATME, DLATMR, DLATMS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, SQRT
      // ..
      // .. Data statements ..
      DATA               KTYPE / 1, 2, 3, 5*4, 4*6, 6*6, 3*9 /
      DATA               KMAGN / 3*1, 1, 1, 1, 2, 3, 4*1, 1, 1, 1, 1, 2, 3, 1, 2, 3 /       DATA               KMODE / 3*0, 4, 3, 1, 4, 4, 4, 3, 1, 5, 4, 3, 1, 5, 5, 5, 4, 3, 1 /
      DATA               KCONDS / 3*0, 5*0, 4*1, 6*2, 3*0 /
      DATA               BAL / 'N', 'P', 'S', 'B' /
      // ..
      // .. Executable Statements ..

      PATH( 1: 1 ) = 'double          ';
      PATH( 2: 3 ) = 'VX'

      // Check for errors

      NTESTT = 0
      NTESTF = 0
      INFO = 0

      // Important constants

      BADNN = .FALSE.

      // 12 is the largest dimension in the input file of precomputed
      // problems

      NMAX = 12
      DO 10 J = 1, NSIZES
         NMAX = MAX( NMAX, NN( J ) )
         IF( NN( J ).LT.0 ) BADNN = .TRUE.
   10 CONTINUE

      // Check for errors

      if ( NSIZES.LT.0 ) {
         INFO = -1
      } else if ( BADNN ) {
         INFO = -2
      } else if ( NTYPES.LT.0 ) {
         INFO = -3
      } else if ( THRESH.LT.ZERO ) {
         INFO = -6
      } else if ( LDA.LT.1 .OR. LDA.LT.NMAX ) {
         INFO = -10
      } else if ( LDVL.LT.1 .OR. LDVL.LT.NMAX ) {
         INFO = -17
      } else if ( LDVR.LT.1 .OR. LDVR.LT.NMAX ) {
         INFO = -19
      } else if ( LDLRE.LT.1 .OR. LDLRE.LT.NMAX ) {
         INFO = -21
      } else if ( 6*NMAX+2*NMAX**2.GT.NWORK ) {
         INFO = -32
      }

      if ( INFO.NE.0 ) {
         xerbla('DDRVVX', -INFO );
         RETURN
      }

      // If nothing to do check on NIUNIT

      IF( NSIZES.EQ.0 .OR. NTYPES.EQ.0 ) GO TO 160

      // More Important constants

      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      ULP = DLAMCH( 'Precision' )
      ULPINV = ONE / ULP
      RTULP = SQRT( ULP )
      RTULPI = ONE / RTULP

      // Loop over sizes, types

      NERRS = 0

      DO 150 JSIZE = 1, NSIZES
         N = NN( JSIZE )
         if ( NSIZES.NE.1 ) {
            MTYPES = MIN( MAXTYP, NTYPES )
         } else {
            MTYPES = MIN( MAXTYP+1, NTYPES )
         }

         DO 140 JTYPE = 1, MTYPES
            IF( .NOT.DOTYPE( JTYPE ) ) GO TO 140

            // Save ISEED in case of an error.

            DO 20 J = 1, 4
               IOLDSD( J ) = ISEED( J )
   20       CONTINUE

            // Compute "A"

            // Control parameters:

            // KMAGN  KCONDS  KMODE        KTYPE
        // =1  O(1)   1       clustered 1  zero
        // =2  large  large   clustered 2  identity
        // =3  small          exponential  Jordan
        // =4                 arithmetic   diagonal, (w/ eigenvalues)
        // =5                 random log   symmetric, w/ eigenvalues
        // =6                 random       general, w/ eigenvalues
        // =7                              random diagonal
        // =8                              random symmetric
        // =9                              random general
        // =10                             random triangular

            IF( MTYPES.GT.MAXTYP ) GO TO 90

            ITYPE = KTYPE( JTYPE )
            IMODE = KMODE( JTYPE )

            // Compute norm

            GO TO ( 30, 40, 50 )KMAGN( JTYPE )

   30       CONTINUE
            ANORM = ONE
            GO TO 60

   40       CONTINUE
            ANORM = OVFL*ULP
            GO TO 60

   50       CONTINUE
            ANORM = UNFL*ULPINV
            GO TO 60

   60       CONTINUE

            dlaset('Full', LDA, N, ZERO, ZERO, A, LDA );
            IINFO = 0
            COND = ULPINV

            // Special Matrices -- Identity & Jordan block

               // Zero

            if ( ITYPE.EQ.1 ) {
               IINFO = 0

            } else if ( ITYPE.EQ.2 ) {

               // Identity

               DO 70 JCOL = 1, N
                  A( JCOL, JCOL ) = ANORM
   70          CONTINUE

            } else if ( ITYPE.EQ.3 ) {

               // Jordan Block

               DO 80 JCOL = 1, N
                  A( JCOL, JCOL ) = ANORM
                  IF( JCOL.GT.1 ) A( JCOL, JCOL-1 ) = ONE
   80          CONTINUE

            } else if ( ITYPE.EQ.4 ) {

               // Diagonal Matrix, [Eigen]values Specified

               dlatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, 0, 0, 'N', A, LDA, WORK( N+1 ), IINFO );

            } else if ( ITYPE.EQ.5 ) {

               // Symmetric, eigenvalues specified

               dlatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, N, N, 'N', A, LDA, WORK( N+1 ), IINFO );

            } else if ( ITYPE.EQ.6 ) {

               // General, eigenvalues specified

               if ( KCONDS( JTYPE ).EQ.1 ) {
                  CONDS = ONE
               } else if ( KCONDS( JTYPE ).EQ.2 ) {
                  CONDS = RTULPI
               } else {
                  CONDS = ZERO
               }

               ADUMMA( 1 ) = ' '
               dlatme(N, 'S', ISEED, WORK, IMODE, COND, ONE, ADUMMA, 'T', 'T', 'T', WORK( N+1 ), 4, CONDS, N, N, ANORM, A, LDA, WORK( 2*N+1 ), IINFO );

            } else if ( ITYPE.EQ.7 ) {

               // Diagonal, random eigenvalues

               dlatmr(N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, 0, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else if ( ITYPE.EQ.8 ) {

               // Symmetric, random eigenvalues

               dlatmr(N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else if ( ITYPE.EQ.9 ) {

               // General, random eigenvalues

               dlatmr(N, N, 'S', ISEED, 'N', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );
               if ( N.GE.4 ) {
                  dlaset('Full', 2, N, ZERO, ZERO, A, LDA );
                  dlaset('Full', N-3, 1, ZERO, ZERO, A( 3, 1 ), LDA )                   CALL DLASET( 'Full', N-3, 2, ZERO, ZERO, A( 3, N-1 ), LDA )                   CALL DLASET( 'Full', 1, N, ZERO, ZERO, A( N, 1 ), LDA );
               }

            } else if ( ITYPE.EQ.10 ) {

               // Triangular, random eigenvalues

               dlatmr(N, N, 'S', ISEED, 'N', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else {

               IINFO = 1
            }

            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9992 )'Generator', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            }

   90       CONTINUE

            // Test for minimal and generous workspace

            DO 130 IWK = 1, 3
               if ( IWK.EQ.1 ) {
                  NNWORK = 3*N
               } else if ( IWK.EQ.2 ) {
                  NNWORK = 6*N + N**2
               } else {
                  NNWORK = 6*N + 2*N**2
               }
               NNWORK = MAX( NNWORK, 1 )

               // Test for all balancing options

               DO 120 IBAL = 1, 4
                  BALANC = BAL( IBAL )

                  // Perform tests

                  dget23(.FALSE., BALANC, JTYPE, THRESH, IOLDSD, NOUNIT, N, A, LDA, H, WR, WI, WR1, WI1, VL, LDVL, VR, LDVR, LRE, LDLRE, RCONDV, RCNDV1, RCDVIN, RCONDE, RCNDE1, RCDEIN, SCALE, SCALE1, RESULT, WORK, NNWORK, IWORK, INFO );

                  // Check for RESULT(j) > THRESH

                  NTEST = 0
                  NFAIL = 0
                  DO 100 J = 1, 9
                     IF( RESULT( J ).GE.ZERO ) NTEST = NTEST + 1                      IF( RESULT( J ).GE.THRESH ) NFAIL = NFAIL + 1
  100             CONTINUE

                  IF( NFAIL.GT.0 ) NTESTF = NTESTF + 1
                  if ( NTESTF.EQ.1 ) {
                     WRITE( NOUNIT, FMT = 9999 )PATH
                     WRITE( NOUNIT, FMT = 9998 )
                     WRITE( NOUNIT, FMT = 9997 )
                     WRITE( NOUNIT, FMT = 9996 )
                     WRITE( NOUNIT, FMT = 9995 )THRESH
                     NTESTF = 2
                  }

                  DO 110 J = 1, 9
                     if ( RESULT( J ).GE.THRESH ) {
                        WRITE( NOUNIT, FMT = 9994 )BALANC, N, IWK, IOLDSD, JTYPE, J, RESULT( J )
                     }
  110             CONTINUE

                  NERRS = NERRS + NFAIL
                  NTESTT = NTESTT + NTEST

  120          CONTINUE
  130       CONTINUE
  140    CONTINUE
  150 CONTINUE

  160 CONTINUE

      // Read in data from file to check accuracy of condition estimation.
      // Assume input eigenvalues are sorted lexicographically (increasing
      // by real part, then decreasing by imaginary part)

      JTYPE = 0
  170 CONTINUE
      READ( NIUNIT, FMT = *, END = 220 )N

      // Read input data until N=0

      IF( N.EQ.0 ) GO TO 220
      JTYPE = JTYPE + 1
      ISEED( 1 ) = JTYPE
      DO 180 I = 1, N
         READ( NIUNIT, FMT = * )( A( I, J ), J = 1, N )
  180 CONTINUE
      DO 190 I = 1, N
         READ( NIUNIT, FMT = * )WR1( I ), WI1( I ), RCDEIN( I ), RCDVIN( I )
  190 CONTINUE
      dget23(.TRUE., 'N', 22, THRESH, ISEED, NOUNIT, N, A, LDA, H, WR, WI, WR1, WI1, VL, LDVL, VR, LDVR, LRE, LDLRE, RCONDV, RCNDV1, RCDVIN, RCONDE, RCNDE1, RCDEIN, SCALE, SCALE1, RESULT, WORK, 6*N+2*N**2, IWORK, INFO );

      // Check for RESULT(j) > THRESH

      NTEST = 0
      NFAIL = 0
      DO 200 J = 1, 11
         IF( RESULT( J ).GE.ZERO ) NTEST = NTEST + 1          IF( RESULT( J ).GE.THRESH ) NFAIL = NFAIL + 1
  200 CONTINUE

      IF( NFAIL.GT.0 ) NTESTF = NTESTF + 1
      if ( NTESTF.EQ.1 ) {
         WRITE( NOUNIT, FMT = 9999 )PATH
         WRITE( NOUNIT, FMT = 9998 )
         WRITE( NOUNIT, FMT = 9997 )
         WRITE( NOUNIT, FMT = 9996 )
         WRITE( NOUNIT, FMT = 9995 )THRESH
         NTESTF = 2
      }

      DO 210 J = 1, 11
         if ( RESULT( J ).GE.THRESH ) {
            WRITE( NOUNIT, FMT = 9993 )N, JTYPE, J, RESULT( J )
         }
  210 CONTINUE

      NERRS = NERRS + NFAIL
      NTESTT = NTESTT + NTEST
      GO TO 170
  220 CONTINUE

      // Summary

      dlasum(PATH, NOUNIT, NERRS, NTESTT );

 9999 FORMAT( / 1X, A3, ' -- Real Eigenvalue-Eigenvector Decomposition', ' Expert Driver', / ' Matrix types (see DDRVVX for details): ' )

 9998 FORMAT( / ' Special Matrices:', / '  1=Zero matrix.             ', '           ', '  5=Diagonal: geometr. spaced entries.', / '  2=Identity matrix.                    ', '  6=Diagona', 'l: clustered entries.', / '  3=Transposed Jordan block.  ', '          ', '  7=Diagonal: large, evenly spaced.', / '  ', '4=Diagonal: evenly spaced entries.    ', '  8=Diagonal: s', 'mall, evenly spaced.' )
 9997 FORMAT( ' Dense, Non-Symmetric Matrices:', / '  9=Well-cond., ev', 'enly spaced eigenvals.', ' 14=Ill-cond., geomet. spaced e', 'igenals.', / ' 10=Well-cond., geom. spaced eigenvals. ', ' 15=Ill-conditioned, clustered e.vals.', / ' 11=Well-cond', 'itioned, clustered e.vals. ', ' 16=Ill-cond., random comp', 'lex ', / ' 12=Well-cond., random complex ', '         ', ' 17=Ill-cond., large rand. complx ', / ' 13=Ill-condi', 'tioned, evenly spaced.     ', ' 18=Ill-cond., small rand.', ' complx ' )
 9996 FORMAT( ' 19=Matrix with random O(1) entries.    ', ' 21=Matrix ', 'with small random entries.', / ' 20=Matrix with large ran', 'dom entries.   ', ' 22=Matrix read from input file', / )
 9995 FORMAT( ' Tests performed with test threshold =', F8.2, / / ' 1 = | A VR - VR W | / ( n |A| ulp ) ', / ' 2 = | transpose(A) VL - VL W | / ( n |A| ulp ) ', / ' 3 = | |VR(i)| - 1 | / ulp ', / ' 4 = | |VL(i)| - 1 | / ulp ', / ' 5 = 0 if W same no matter if VR or VL computed,', ' 1/ulp otherwise', / ' 6 = 0 if VR same no matter what else computed,', '  1/ulp otherwise', / ' 7 = 0 if VL same no matter what else computed,', '  1/ulp otherwise', / ' 8 = 0 if RCONDV same no matter what else computed,', '  1/ulp otherwise', / ' 9 = 0 if SCALE, ILO, IHI, ABNRM same no matter what else', ' computed,  1/ulp otherwise', / ' 10 = | RCONDV - RCONDV(precomputed) | / cond(RCONDV),', / ' 11 = | RCONDE - RCONDE(precomputed) | / cond(RCONDE),' )
 9994 FORMAT( ' BALANC=''', A1, ''',N=', I4, ',IWK=', I1, ', seed=', 4( I4, ',' ), ' type ', I2, ', test(', I2, ')=', G10.3 )
 9993 FORMAT( ' N=', I5, ', input example =', I3, ',  test(', I2, ')=', G10.3 )
 9992 FORMAT( ' DDRVVX: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )

      RETURN

      // End of DDRVVX

      }
