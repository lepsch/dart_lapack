      SUBROUTINE DCHKHS( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH, NOUNIT, A, LDA, H, T1, T2, U, LDU, Z, UZ, WR1, WI1, WR2, WI2, WR3, WI3, EVECTL, EVECTR, EVECTY, EVECTX, UU, TAU, WORK, NWORK, IWORK, SELECT, RESULT, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDU, NOUNIT, NSIZES, NTYPES, NWORK;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * ), SELECT( * );
      int                ISEED( 4 ), IWORK( * ), NN( * );
      double             A( LDA, * ), EVECTL( LDU, * ), EVECTR( LDU, * ), EVECTX( LDU, * ), EVECTY( LDU, * ), H( LDA, * ), RESULT( 16 ), T1( LDA, * ), T2( LDA, * ), TAU( * ), U( LDU, * ), UU( LDU, * ), UZ( LDU, * ), WI1( * ), WI2( * ), WI3( * ), WORK( * ), WR1( * ), WR2( * ), WR3( * ), Z( LDU, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
      int                MAXTYP;
      const              MAXTYP = 21 ;
      // ..
      // .. Local Scalars ..
      bool               BADNN, MATCH;
      int                I, IHI, IINFO, ILO, IMODE, IN, ITYPE, J, JCOL, JJ, JSIZE, JTYPE, K, MTYPES, N, N1, NERRS, NMATS, NMAX, NSELC, NSELR, NTEST, NTESTT;
      double             ANINV, ANORM, COND, CONDS, OVFL, RTOVFL, RTULP, RTULPI, RTUNFL, TEMP1, TEMP2, ULP, ULPINV, UNFL;
      // ..
      // .. Local Arrays ..
      String             ADUMMA( 1 );
      int                IDUMMA( 1 ), IOLDSD( 4 ), KCONDS( MAXTYP ), KMAGN( MAXTYP ), KMODE( MAXTYP ), KTYPE( MAXTYP );
      double             DUMMA( 6 );
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DGEHRD, DGEMM, DGET10, DGET22, DHSEIN, DHSEQR, DHST01, DLACPY, DLAFTS, DLASET, DLASUM, DLATME, DLATMR, DLATMS, DORGHR, DORMHR, DTREVC, DTREVC3, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, MIN, SQRT
      // ..
      // .. Data statements ..
      DATA               KTYPE / 1, 2, 3, 5*4, 4*6, 6*6, 3*9 /
      DATA               KMAGN / 3*1, 1, 1, 1, 2, 3, 4*1, 1, 1, 1, 1, 2, 3, 1, 2, 3 /       DATA               KMODE / 3*0, 4, 3, 1, 4, 4, 4, 3, 1, 5, 4, 3, 1, 5, 5, 5, 4, 3, 1 /
      DATA               KCONDS / 3*0, 5*0, 4*1, 6*2, 3*0 /
      // ..
      // .. Executable Statements ..

      // Check for errors

      NTESTT = 0
      INFO = 0

      BADNN = .FALSE.
      NMAX = 0
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
      } else if ( LDA.LE.1 .OR. LDA.LT.NMAX ) {
         INFO = -9
      } else if ( LDU.LE.1 .OR. LDU.LT.NMAX ) {
         INFO = -14
      } else if ( 4*NMAX*NMAX+2.GT.NWORK ) {
         INFO = -28
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'DCHKHS', -INFO )
         RETURN
      }

      // Quick return if possible

      IF( NSIZES.EQ.0 .OR. NTYPES.EQ.0 ) RETURN

      // More important constants

      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = DLAMCH( 'Overflow' )
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )
      ULPINV = ONE / ULP
      RTUNFL = SQRT( UNFL )
      RTOVFL = SQRT( OVFL )
      RTULP = SQRT( ULP )
      RTULPI = ONE / RTULP

      // Loop over sizes, types

      NERRS = 0
      NMATS = 0

      DO 270 JSIZE = 1, NSIZES
         N = NN( JSIZE )
         IF( N.EQ.0 ) GO TO 270
         N1 = MAX( 1, N )
         ANINV = ONE / DBLE( N1 )

         if ( NSIZES.NE.1 ) {
            MTYPES = MIN( MAXTYP, NTYPES )
         } else {
            MTYPES = MIN( MAXTYP+1, NTYPES )
         }

         DO 260 JTYPE = 1, MTYPES
            IF( .NOT.DOTYPE( JTYPE ) ) GO TO 260
            NMATS = NMATS + 1
            NTEST = 0

            // Save ISEED in case of an error.

            DO 20 J = 1, 4
               IOLDSD( J ) = ISEED( J )
   20       CONTINUE

            // Initialize RESULT

            DO 30 J = 1, 16
               RESULT( J ) = ZERO
   30       CONTINUE

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

            IF( MTYPES.GT.MAXTYP ) GO TO 100

            ITYPE = KTYPE( JTYPE )
            IMODE = KMODE( JTYPE )

            // Compute norm

            GO TO ( 40, 50, 60 )KMAGN( JTYPE )

   40       CONTINUE
            ANORM = ONE
            GO TO 70

   50       CONTINUE
            ANORM = ( RTOVFL*ULP )*ANINV
            GO TO 70

   60       CONTINUE
            ANORM = RTUNFL*N*ULPINV
            GO TO 70

   70       CONTINUE

            CALL DLASET( 'Full', LDA, N, ZERO, ZERO, A, LDA )
            IINFO = 0
            COND = ULPINV

            // Special Matrices

            if ( ITYPE.EQ.1 ) {

               // Zero

               IINFO = 0

            } else if ( ITYPE.EQ.2 ) {

               // Identity

               DO 80 JCOL = 1, N
                  A( JCOL, JCOL ) = ANORM
   80          CONTINUE

            } else if ( ITYPE.EQ.3 ) {

               // Jordan Block

               DO 90 JCOL = 1, N
                  A( JCOL, JCOL ) = ANORM
                  IF( JCOL.GT.1 ) A( JCOL, JCOL-1 ) = ONE
   90          CONTINUE

            } else if ( ITYPE.EQ.4 ) {

               // Diagonal Matrix, [Eigen]values Specified

               CALL DLATMS( N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, 0, 0, 'N', A, LDA, WORK( N+1 ), IINFO )

            } else if ( ITYPE.EQ.5 ) {

               // Symmetric, eigenvalues specified

               CALL DLATMS( N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, N, N, 'N', A, LDA, WORK( N+1 ), IINFO )

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
               CALL DLATME( N, 'S', ISEED, WORK, IMODE, COND, ONE, ADUMMA, 'T', 'T', 'T', WORK( N+1 ), 4, CONDS, N, N, ANORM, A, LDA, WORK( 2*N+1 ), IINFO )

            } else if ( ITYPE.EQ.7 ) {

               // Diagonal, random eigenvalues

               CALL DLATMR( N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, 0, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )

            } else if ( ITYPE.EQ.8 ) {

               // Symmetric, random eigenvalues

               CALL DLATMR( N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )

            } else if ( ITYPE.EQ.9 ) {

               // General, random eigenvalues

               CALL DLATMR( N, N, 'S', ISEED, 'N', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )

            } else if ( ITYPE.EQ.10 ) {

               // Triangular, random eigenvalues

               CALL DLATMR( N, N, 'S', ISEED, 'N', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )

            } else {

               IINFO = 1
            }

            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'Generator', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            }

  100       CONTINUE

            // Call DGEHRD to compute H and U, do tests.

            CALL DLACPY( ' ', N, N, A, LDA, H, LDA )

            NTEST = 1

            ILO = 1
            IHI = N

            CALL DGEHRD( N, ILO, IHI, H, LDA, WORK, WORK( N+1 ), NWORK-N, IINFO )

            if ( IINFO.NE.0 ) {
               RESULT( 1 ) = ULPINV
               WRITE( NOUNIT, FMT = 9999 )'DGEHRD', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               GO TO 250
            }

            DO 120 J = 1, N - 1
               UU( J+1, J ) = ZERO
               DO 110 I = J + 2, N
                  U( I, J ) = H( I, J )
                  UU( I, J ) = H( I, J )
                  H( I, J ) = ZERO
  110          CONTINUE
  120       CONTINUE
            CALL DCOPY( N-1, WORK, 1, TAU, 1 )
            CALL DORGHR( N, ILO, IHI, U, LDU, WORK, WORK( N+1 ), NWORK-N, IINFO )
            NTEST = 2

            CALL DHST01( N, ILO, IHI, A, LDA, H, LDA, U, LDU, WORK, NWORK, RESULT( 1 ) )

            // Call DHSEQR to compute T1, T2 and Z, do tests.

            // Eigenvalues only (WR3,WI3)

            CALL DLACPY( ' ', N, N, H, LDA, T2, LDA )
            NTEST = 3
            RESULT( 3 ) = ULPINV

            CALL DHSEQR( 'E', 'N', N, ILO, IHI, T2, LDA, WR3, WI3, UZ, LDU, WORK, NWORK, IINFO )
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'DHSEQR(E)', IINFO, N, JTYPE, IOLDSD
               if ( IINFO.LE.N+2 ) {
                  INFO = ABS( IINFO )
                  GO TO 250
               }
            }

            // Eigenvalues (WR2,WI2) and Full Schur Form (T2)

            CALL DLACPY( ' ', N, N, H, LDA, T2, LDA )

            CALL DHSEQR( 'S', 'N', N, ILO, IHI, T2, LDA, WR2, WI2, UZ, LDU, WORK, NWORK, IINFO )
            if ( IINFO.NE.0 .AND. IINFO.LE.N+2 ) {
               WRITE( NOUNIT, FMT = 9999 )'DHSEQR(S)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               GO TO 250
            }

            // Eigenvalues (WR1,WI1), Schur Form (T1), and Schur vectors
            // (UZ)

            CALL DLACPY( ' ', N, N, H, LDA, T1, LDA )
            CALL DLACPY( ' ', N, N, U, LDU, UZ, LDU )

            CALL DHSEQR( 'S', 'V', N, ILO, IHI, T1, LDA, WR1, WI1, UZ, LDU, WORK, NWORK, IINFO )
            if ( IINFO.NE.0 .AND. IINFO.LE.N+2 ) {
               WRITE( NOUNIT, FMT = 9999 )'DHSEQR(V)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               GO TO 250
            }

            // Compute Z = U' UZ

            CALL DGEMM( 'T', 'N', N, N, N, ONE, U, LDU, UZ, LDU, ZERO, Z, LDU )
            NTEST = 8

            // Do Tests 3: | H - Z T Z' | / ( |H| n ulp )
                 // and 4: | I - Z Z' | / ( n ulp )

            CALL DHST01( N, ILO, IHI, H, LDA, T1, LDA, Z, LDU, WORK, NWORK, RESULT( 3 ) )

            // Do Tests 5: | A - UZ T (UZ)' | / ( |A| n ulp )
                 // and 6: | I - UZ (UZ)' | / ( n ulp )

            CALL DHST01( N, ILO, IHI, A, LDA, T1, LDA, UZ, LDU, WORK, NWORK, RESULT( 5 ) )

            // Do Test 7: | T2 - T1 | / ( |T| n ulp )

            CALL DGET10( N, N, T2, LDA, T1, LDA, WORK, RESULT( 7 ) )

            // Do Test 8: | W2 - W1 | / ( max(|W1|,|W2|) ulp )

            TEMP1 = ZERO
            TEMP2 = ZERO
            DO 130 J = 1, N
               TEMP1 = MAX( TEMP1, ABS( WR1( J ) )+ABS( WI1( J ) ), ABS( WR2( J ) )+ABS( WI2( J ) ) )                TEMP2 = MAX( TEMP2, ABS( WR1( J )-WR2( J ) )+ ABS( WI1( J )-WI2( J ) ) )
  130       CONTINUE

            RESULT( 8 ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )

            // Compute the Left and Right Eigenvectors of T

            // Compute the Right eigenvector Matrix:

            NTEST = 9
            RESULT( 9 ) = ULPINV

            // Select last max(N/4,1) real, max(N/4,1) complex eigenvectors

            NSELC = 0
            NSELR = 0
            J = N
  140       CONTINUE
            if ( WI1( J ).EQ.ZERO ) {
               if ( NSELR.LT.MAX( N / 4, 1 ) ) {
                  NSELR = NSELR + 1
                  SELECT( J ) = .TRUE.
               } else {
                  SELECT( J ) = .FALSE.
               }
               J = J - 1
            } else {
               if ( NSELC.LT.MAX( N / 4, 1 ) ) {
                  NSELC = NSELC + 1
                  SELECT( J ) = .TRUE.
                  SELECT( J-1 ) = .FALSE.
               } else {
                  SELECT( J ) = .FALSE.
                  SELECT( J-1 ) = .FALSE.
               }
               J = J - 2
            }
            IF( J.GT.0 ) GO TO 140

            CALL DTREVC( 'Right', 'All', SELECT, N, T1, LDA, DUMMA, LDU, EVECTR, LDU, N, IN, WORK, IINFO )
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'DTREVC(R,A)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               GO TO 250
            }

            // Test 9:  | TR - RW | / ( |T| |R| ulp )

            CALL DGET22( 'N', 'N', 'N', N, T1, LDA, EVECTR, LDU, WR1, WI1, WORK, DUMMA( 1 ) )
            RESULT( 9 ) = DUMMA( 1 )
            if ( DUMMA( 2 ).GT.THRESH ) {
               WRITE( NOUNIT, FMT = 9998 )'Right', 'DTREVC', DUMMA( 2 ), N, JTYPE, IOLDSD
            }

            // Compute selected right eigenvectors and confirm that
            // they agree with previous right eigenvectors

            CALL DTREVC( 'Right', 'Some', SELECT, N, T1, LDA, DUMMA, LDU, EVECTL, LDU, N, IN, WORK, IINFO )
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'DTREVC(R,S)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               GO TO 250
            }

            K = 1
            MATCH = .TRUE.
            DO 170 J = 1, N
               if ( SELECT( J ) .AND. WI1( J ).EQ.ZERO ) {
                  DO 150 JJ = 1, N
                     if ( EVECTR( JJ, J ).NE.EVECTL( JJ, K ) ) {
                        MATCH = .FALSE.
                        GO TO 180
                     }
  150             CONTINUE
                  K = K + 1
               } else if ( SELECT( J ) .AND. WI1( J ).NE.ZERO ) {
                  DO 160 JJ = 1, N
                     if ( EVECTR( JJ, J ).NE.EVECTL( JJ, K ) .OR. EVECTR( JJ, J+1 ).NE.EVECTL( JJ, K+1 ) ) {
                        MATCH = .FALSE.
                        GO TO 180
                     }
  160             CONTINUE
                  K = K + 2
               }
  170       CONTINUE
  180       CONTINUE
            IF( .NOT.MATCH ) WRITE( NOUNIT, FMT = 9997 )'Right', 'DTREVC', N, JTYPE, IOLDSD

            // Compute the Left eigenvector Matrix:

            NTEST = 10
            RESULT( 10 ) = ULPINV
            CALL DTREVC( 'Left', 'All', SELECT, N, T1, LDA, EVECTL, LDU, DUMMA, LDU, N, IN, WORK, IINFO )
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'DTREVC(L,A)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               GO TO 250
            }

            // Test 10:  | LT - WL | / ( |T| |L| ulp )

            CALL DGET22( 'Trans', 'N', 'Conj', N, T1, LDA, EVECTL, LDU, WR1, WI1, WORK, DUMMA( 3 ) )
            RESULT( 10 ) = DUMMA( 3 )
            if ( DUMMA( 4 ).GT.THRESH ) {
               WRITE( NOUNIT, FMT = 9998 )'Left', 'DTREVC', DUMMA( 4 ), N, JTYPE, IOLDSD
            }

            // Compute selected left eigenvectors and confirm that
            // they agree with previous left eigenvectors

            CALL DTREVC( 'Left', 'Some', SELECT, N, T1, LDA, EVECTR, LDU, DUMMA, LDU, N, IN, WORK, IINFO )
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'DTREVC(L,S)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               GO TO 250
            }

            K = 1
            MATCH = .TRUE.
            DO 210 J = 1, N
               if ( SELECT( J ) .AND. WI1( J ).EQ.ZERO ) {
                  DO 190 JJ = 1, N
                     if ( EVECTL( JJ, J ).NE.EVECTR( JJ, K ) ) {
                        MATCH = .FALSE.
                        GO TO 220
                     }
  190             CONTINUE
                  K = K + 1
               } else if ( SELECT( J ) .AND. WI1( J ).NE.ZERO ) {
                  DO 200 JJ = 1, N
                     if ( EVECTL( JJ, J ).NE.EVECTR( JJ, K ) .OR. EVECTL( JJ, J+1 ).NE.EVECTR( JJ, K+1 ) ) {
                        MATCH = .FALSE.
                        GO TO 220
                     }
  200             CONTINUE
                  K = K + 2
               }
  210       CONTINUE
  220       CONTINUE
            IF( .NOT.MATCH ) WRITE( NOUNIT, FMT = 9997 )'Left', 'DTREVC', N, JTYPE, IOLDSD

            // Call DHSEIN for Right eigenvectors of H, do test 11

            NTEST = 11
            RESULT( 11 ) = ULPINV
            DO 230 J = 1, N
               SELECT( J ) = .TRUE.
  230       CONTINUE

            CALL DHSEIN( 'Right', 'Qr', 'Ninitv', SELECT, N, H, LDA, WR3, WI3, DUMMA, LDU, EVECTX, LDU, N1, IN, WORK, IWORK, IWORK, IINFO )
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'DHSEIN(R)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO.LT.0 ) GO TO 250
            } else {

               // Test 11:  | HX - XW | / ( |H| |X| ulp )

                         // (from inverse iteration)

               CALL DGET22( 'N', 'N', 'N', N, H, LDA, EVECTX, LDU, WR3, WI3, WORK, DUMMA( 1 ) )                IF( DUMMA( 1 ).LT.ULPINV ) RESULT( 11 ) = DUMMA( 1 )*ANINV
               if ( DUMMA( 2 ).GT.THRESH ) {
                  WRITE( NOUNIT, FMT = 9998 )'Right', 'DHSEIN', DUMMA( 2 ), N, JTYPE, IOLDSD
               }
            }

            // Call DHSEIN for Left eigenvectors of H, do test 12

            NTEST = 12
            RESULT( 12 ) = ULPINV
            DO 240 J = 1, N
               SELECT( J ) = .TRUE.
  240       CONTINUE

            CALL DHSEIN( 'Left', 'Qr', 'Ninitv', SELECT, N, H, LDA, WR3, WI3, EVECTY, LDU, DUMMA, LDU, N1, IN, WORK, IWORK, IWORK, IINFO )
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'DHSEIN(L)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO.LT.0 ) GO TO 250
            } else {

               // Test 12:  | YH - WY | / ( |H| |Y| ulp )

                         // (from inverse iteration)

               CALL DGET22( 'C', 'N', 'C', N, H, LDA, EVECTY, LDU, WR3, WI3, WORK, DUMMA( 3 ) )                IF( DUMMA( 3 ).LT.ULPINV ) RESULT( 12 ) = DUMMA( 3 )*ANINV
               if ( DUMMA( 4 ).GT.THRESH ) {
                  WRITE( NOUNIT, FMT = 9998 )'Left', 'DHSEIN', DUMMA( 4 ), N, JTYPE, IOLDSD
               }
            }

            // Call DORMHR for Right eigenvectors of A, do test 13

            NTEST = 13
            RESULT( 13 ) = ULPINV

            CALL DORMHR( 'Left', 'No transpose', N, N, ILO, IHI, UU, LDU, TAU, EVECTX, LDU, WORK, NWORK, IINFO )
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'DORMHR(R)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO.LT.0 ) GO TO 250
            } else {

               // Test 13:  | AX - XW | / ( |A| |X| ulp )

                         // (from inverse iteration)

               CALL DGET22( 'N', 'N', 'N', N, A, LDA, EVECTX, LDU, WR3, WI3, WORK, DUMMA( 1 ) )                IF( DUMMA( 1 ).LT.ULPINV ) RESULT( 13 ) = DUMMA( 1 )*ANINV
            }

            // Call DORMHR for Left eigenvectors of A, do test 14

            NTEST = 14
            RESULT( 14 ) = ULPINV

            CALL DORMHR( 'Left', 'No transpose', N, N, ILO, IHI, UU, LDU, TAU, EVECTY, LDU, WORK, NWORK, IINFO )
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'DORMHR(L)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               IF( IINFO.LT.0 ) GO TO 250
            } else {

               // Test 14:  | YA - WY | / ( |A| |Y| ulp )

                         // (from inverse iteration)

               CALL DGET22( 'C', 'N', 'C', N, A, LDA, EVECTY, LDU, WR3, WI3, WORK, DUMMA( 3 ) )                IF( DUMMA( 3 ).LT.ULPINV ) RESULT( 14 ) = DUMMA( 3 )*ANINV
            }

            // Compute Left and Right Eigenvectors of A

            // Compute a Right eigenvector matrix:

            NTEST = 15
            RESULT( 15 ) = ULPINV

            CALL DLACPY( ' ', N, N, UZ, LDU, EVECTR, LDU )

            CALL DTREVC3( 'Right', 'Back', SELECT, N, T1, LDA, DUMMA, LDU, EVECTR, LDU, N, IN, WORK, NWORK, IINFO )
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'DTREVC3(R,B)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               GO TO 250
            }

            // Test 15:  | AR - RW | / ( |A| |R| ulp )

                      // (from Schur decomposition)

            CALL DGET22( 'N', 'N', 'N', N, A, LDA, EVECTR, LDU, WR1, WI1, WORK, DUMMA( 1 ) )
            RESULT( 15 ) = DUMMA( 1 )
            if ( DUMMA( 2 ).GT.THRESH ) {
               WRITE( NOUNIT, FMT = 9998 )'Right', 'DTREVC3', DUMMA( 2 ), N, JTYPE, IOLDSD
            }

            // Compute a Left eigenvector matrix:

            NTEST = 16
            RESULT( 16 ) = ULPINV

            CALL DLACPY( ' ', N, N, UZ, LDU, EVECTL, LDU )

            CALL DTREVC3( 'Left', 'Back', SELECT, N, T1, LDA, EVECTL, LDU, DUMMA, LDU, N, IN, WORK, NWORK, IINFO )
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'DTREVC3(L,B)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               GO TO 250
            }

            // Test 16:  | LA - WL | / ( |A| |L| ulp )

                      // (from Schur decomposition)

            CALL DGET22( 'Trans', 'N', 'Conj', N, A, LDA, EVECTL, LDU, WR1, WI1, WORK, DUMMA( 3 ) )
            RESULT( 16 ) = DUMMA( 3 )
            if ( DUMMA( 4 ).GT.THRESH ) {
               WRITE( NOUNIT, FMT = 9998 )'Left', 'DTREVC3', DUMMA( 4 ), N, JTYPE, IOLDSD
            }

            // End of Loop -- Check for RESULT(j) > THRESH

  250       CONTINUE

            NTESTT = NTESTT + NTEST
            CALL DLAFTS( 'DHS', N, N, JTYPE, NTEST, RESULT, IOLDSD, THRESH, NOUNIT, NERRS )

  260    CONTINUE
  270 CONTINUE

      // Summary

      CALL DLASUM( 'DHS', NOUNIT, NERRS, NTESTT )

      RETURN

 9999 FORMAT( ' DCHKHS: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )
 9998 FORMAT( ' DCHKHS: ', A, ' Eigenvectors from ', A, ' incorrectly ', 'normalized.', / ' Bits of error=', 0P, G10.3, ',', 9X, 'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )
 9997 FORMAT( ' DCHKHS: Selected ', A, ' Eigenvectors from ', A, ' do not match other eigenvectors ', 9X, 'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )

      // End of DCHKHS

      }
