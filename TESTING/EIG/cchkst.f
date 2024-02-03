      SUBROUTINE CCHKST( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH, NOUNIT, A, LDA, AP, SD, SE, D1, D2, D3, D4, D5, WA1, WA2, WA3, WR, U, LDU, V, VP, TAU, Z, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, RESULT, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDU, LIWORK, LRWORK, LWORK, NOUNIT, NSIZES, NTYPES;
      REAL               THRESH
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                ISEED( 4 ), IWORK( * ), NN( * );
      REAL               D1( * ), D2( * ), D3( * ), D4( * ), D5( * ), RESULT( * ), RWORK( * ), SD( * ), SE( * ), WA1( * ), WA2( * ), WA3( * ), WR( * )
      COMPLEX            A( LDA, * ), AP( * ), TAU( * ), U( LDU, * ), V( LDU, * ), VP( * ), WORK( * ), Z( LDU, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TWO, EIGHT, TEN, HUN
      const              ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0, EIGHT = 8.0E0, TEN = 10.0E0, HUN = 100.0E0 ;
      COMPLEX            CZERO, CONE
      const              CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) ;
      REAL               HALF
      const              HALF = ONE / TWO ;
      int                MAXTYP;
      const              MAXTYP = 21 ;
      bool               CRANGE;
      const              CRANGE = false ;
      bool               CREL;
      const              CREL = false ;
      // ..
      // .. Local Scalars ..
      bool               BADNN, TRYRAC;
      int                I, IINFO, IL, IMODE, INDE, INDRWK, ITEMP, ITYPE, IU, J, JC, JR, JSIZE, JTYPE, LGN, LIWEDC, LOG2UI, LRWEDC, LWEDC, M, M2, M3, MTYPES, N, NAP, NBLOCK, NERRS, NMATS, NMAX, NSPLIT, NTEST, NTESTT;
      REAL               ABSTOL, ANINV, ANORM, COND, OVFL, RTOVFL, RTUNFL, TEMP1, TEMP2, TEMP3, TEMP4, ULP, ULPINV, UNFL, VL, VU
      // ..
      // .. Local Arrays ..
      int                IDUMMA( 1 ), IOLDSD( 4 ), ISEED2( 4 ), KMAGN( MAXTYP ), KMODE( MAXTYP ), KTYPE( MAXTYP );
      REAL               DUMMA( 1 )
      // ..
      // .. External Functions ..
      int                ILAENV;
      REAL               SLAMCH, SLARND, SSXT1
      // EXTERNAL ILAENV, SLAMCH, SLARND, SSXT1
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CHET21, CHETRD, CHPT21, CHPTRD, CLACPY, CLASET, CLATMR, CLATMS, CPTEQR, CSTEDC, CSTEMR, CSTEIN, CSTEQR, CSTT21, CSTT22, CUNGTR, CUPGTR, SCOPY, SLASUM, SSTEBZ, SSTECH, SSTERF, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, CONJG, INT, LOG, MAX, MIN, REAL, SQRT
      // ..
      // .. Data statements ..
      DATA               KTYPE / 1, 2, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 8, 8, 8, 9, 9, 9, 9, 9, 10 /
      DATA               KMAGN / 1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3, 1, 1, 1, 2, 3, 1 /
      DATA               KMODE / 0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0, 4, 3, 1, 4, 4, 3 /
      // ..
      // .. Executable Statements ..

      // Keep ftnchek happy
      IDUMMA( 1 ) = 1

      // Check for errors

      NTESTT = 0
      INFO = 0

      // Important constants

      BADNN = false;
      TRYRAC = true;
      NMAX = 1
      for (J = 1; J <= NSIZES; J++) { // 10
         NMAX = MAX( NMAX, NN( J ) )
         IF( NN( J ).LT.0 ) BADNN = true;
      } // 10

      NBLOCK = ILAENV( 1, 'CHETRD', 'L', NMAX, -1, -1, -1 )
      NBLOCK = MIN( NMAX, MAX( 1, NBLOCK ) )

      // Check for errors

      if ( NSIZES.LT.0 ) {
         INFO = -1
      } else if ( BADNN ) {
         INFO = -2
      } else if ( NTYPES.LT.0 ) {
         INFO = -3
      } else if ( LDA.LT.NMAX ) {
         INFO = -9
      } else if ( LDU.LT.NMAX ) {
         INFO = -23
      } else if ( 2*MAX( 2, NMAX )**2.GT.LWORK ) {
         INFO = -29
      }

      if ( INFO.NE.0 ) {
         xerbla('CCHKST', -INFO );
         RETURN
      }

      // Quick return if possible

      if (NSIZES == 0 .OR. NTYPES == 0) RETURN;

      // More Important constants

      UNFL = SLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )
      ULPINV = ONE / ULP
      LOG2UI = INT( LOG( ULPINV ) / LOG( TWO ) )
      RTUNFL = SQRT( UNFL )
      RTOVFL = SQRT( OVFL )

      // Loop over sizes, types

      for (I = 1; I <= 4; I++) { // 20
         ISEED2( I ) = ISEED( I )
      } // 20
      NERRS = 0
      NMATS = 0

      for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) { // 310
         N = NN( JSIZE )
         if ( N.GT.0 ) {
            LGN = INT( LOG( REAL( N ) ) / LOG( TWO ) )
            if (2**LGN.LT.N) LGN = LGN + 1             IF( 2**LGN.LT.N ) LGN = LGN + 1;
            LWEDC = 1 + 4*N + 2*N*LGN + 4*N**2
            LRWEDC = 1 + 3*N + 2*N*LGN + 4*N**2
            LIWEDC = 6 + 6*N + 5*N*LGN
         } else {
            LWEDC = 8
            LRWEDC = 7
            LIWEDC = 12
         }
         NAP = ( N*( N+1 ) ) / 2
         ANINV = ONE / REAL( MAX( 1, N ) )

         if ( NSIZES.NE.1 ) {
            MTYPES = MIN( MAXTYP, NTYPES )
         } else {
            MTYPES = MIN( MAXTYP+1, NTYPES )
         }

         for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) { // 300
            IF( .NOT.DOTYPE( JTYPE ) ) GO TO 300
            NMATS = NMATS + 1
            NTEST = 0

            for (J = 1; J <= 4; J++) { // 30
               IOLDSD( J ) = ISEED( J )
            } // 30

            // Compute "A"

            // Control parameters:

                // KMAGN  KMODE        KTYPE
            // =1  O(1)   clustered 1  zero
            // =2  large  clustered 2  identity
            // =3  small  exponential  (none)
            // =4         arithmetic   diagonal, (w/ eigenvalues)
            // =5         random log   Hermitian, w/ eigenvalues
            // =6         random       (none)
            // =7                      random diagonal
            // =8                      random Hermitian
            // =9                      positive definite
            // =10                     diagonally dominant tridiagonal

            if (MTYPES.GT.MAXTYP) GO TO 100;

            ITYPE = KTYPE( JTYPE )
            IMODE = KMODE( JTYPE )

            // Compute norm

            GO TO ( 40, 50, 60 )KMAGN( JTYPE )

            } // 40
            ANORM = ONE
            GO TO 70

            } // 50
            ANORM = ( RTOVFL*ULP )*ANINV
            GO TO 70

            } // 60
            ANORM = RTUNFL*N*ULPINV
            GO TO 70

            } // 70

            claset('Full', LDA, N, CZERO, CZERO, A, LDA );
            IINFO = 0
            if ( JTYPE.LE.15 ) {
               COND = ULPINV
            } else {
               COND = ULPINV*ANINV / TEN
            }

            // Special Matrices -- Identity & Jordan block

               // Zero

            if ( ITYPE == 1 ) {
               IINFO = 0

            } else if ( ITYPE == 2 ) {

               // Identity

               for (JC = 1; JC <= N; JC++) { // 80
                  A( JC, JC ) = ANORM
               } // 80

            } else if ( ITYPE == 4 ) {

               // Diagonal Matrix, [Eigen]values Specified

               clatms(N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, 0, 0, 'N', A, LDA, WORK, IINFO );


            } else if ( ITYPE == 5 ) {

               // Hermitian, eigenvalues specified

               clatms(N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, N, N, 'N', A, LDA, WORK, IINFO );

            } else if ( ITYPE == 7 ) {

               // Diagonal, random eigenvalues

               clatmr(N, N, 'S', ISEED, 'H', WORK, 6, ONE, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, 0, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else if ( ITYPE == 8 ) {

               // Hermitian, random eigenvalues

               clatmr(N, N, 'S', ISEED, 'H', WORK, 6, ONE, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else if ( ITYPE == 9 ) {

               // Positive definite, eigenvalues specified.

               clatms(N, N, 'S', ISEED, 'P', RWORK, IMODE, COND, ANORM, N, N, 'N', A, LDA, WORK, IINFO );

            } else if ( ITYPE == 10 ) {

               // Positive definite tridiagonal, eigenvalues specified.

               clatms(N, N, 'S', ISEED, 'P', RWORK, IMODE, COND, ANORM, 1, 1, 'N', A, LDA, WORK, IINFO );
               for (I = 2; I <= N; I++) { // 90
                  TEMP1 = ABS( A( I-1, I ) )
                  TEMP2 = SQRT( ABS( A( I-1, I-1 )*A( I, I ) ) )
                  if ( TEMP1.GT.HALF*TEMP2 ) {
                     A( I-1, I ) = A( I-1, I )* ( HALF*TEMP2 / ( UNFL+TEMP1 ) )
                     A( I, I-1 ) = CONJG( A( I-1, I ) )
                  }
               } // 90

            } else {

               IINFO = 1
            }

            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'Generator', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            }

            } // 100

            // Call CHETRD and CUNGTR to compute S and U from
            // upper triangle.

            clacpy('U', N, N, A, LDA, V, LDU );

            NTEST = 1
            chetrd('U', N, V, LDU, SD, SE, TAU, WORK, LWORK, IINFO );

            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'CHETRD(U)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 1 ) = ULPINV
                  GO TO 280
               }
            }

            clacpy('U', N, N, V, LDU, U, LDU );

            NTEST = 2
            cungtr('U', N, U, LDU, TAU, WORK, LWORK, IINFO );
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'CUNGTR(U)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 2 ) = ULPINV
                  GO TO 280
               }
            }

            // Do tests 1 and 2

            chet21(2, 'Upper', N, 1, A, LDA, SD, SE, U, LDU, V, LDU, TAU, WORK, RWORK, RESULT( 1 ) );
            chet21(3, 'Upper', N, 1, A, LDA, SD, SE, U, LDU, V, LDU, TAU, WORK, RWORK, RESULT( 2 ) );

            // Call CHETRD and CUNGTR to compute S and U from
            // lower triangle, do tests.

            clacpy('L', N, N, A, LDA, V, LDU );

            NTEST = 3
            chetrd('L', N, V, LDU, SD, SE, TAU, WORK, LWORK, IINFO );

            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'CHETRD(L)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 3 ) = ULPINV
                  GO TO 280
               }
            }

            clacpy('L', N, N, V, LDU, U, LDU );

            NTEST = 4
            cungtr('L', N, U, LDU, TAU, WORK, LWORK, IINFO );
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'CUNGTR(L)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 4 ) = ULPINV
                  GO TO 280
               }
            }

            chet21(2, 'Lower', N, 1, A, LDA, SD, SE, U, LDU, V, LDU, TAU, WORK, RWORK, RESULT( 3 ) );
            chet21(3, 'Lower', N, 1, A, LDA, SD, SE, U, LDU, V, LDU, TAU, WORK, RWORK, RESULT( 4 ) );

            // Store the upper triangle of A in AP

            I = 0
            for (JC = 1; JC <= N; JC++) { // 120
               for (JR = 1; JR <= JC; JR++) { // 110
                  I = I + 1
                  AP( I ) = A( JR, JC )
               } // 110
            } // 120

            // Call CHPTRD and CUPGTR to compute S and U from AP

            ccopy(NAP, AP, 1, VP, 1 );

            NTEST = 5
            chptrd('U', N, VP, SD, SE, TAU, IINFO );

            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'CHPTRD(U)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 5 ) = ULPINV
                  GO TO 280
               }
            }

            NTEST = 6
            cupgtr('U', N, VP, TAU, U, LDU, WORK, IINFO );
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'CUPGTR(U)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 6 ) = ULPINV
                  GO TO 280
               }
            }

            // Do tests 5 and 6

            chpt21(2, 'Upper', N, 1, AP, SD, SE, U, LDU, VP, TAU, WORK, RWORK, RESULT( 5 ) );
            chpt21(3, 'Upper', N, 1, AP, SD, SE, U, LDU, VP, TAU, WORK, RWORK, RESULT( 6 ) );

            // Store the lower triangle of A in AP

            I = 0
            for (JC = 1; JC <= N; JC++) { // 140
               for (JR = JC; JR <= N; JR++) { // 130
                  I = I + 1
                  AP( I ) = A( JR, JC )
               } // 130
            } // 140

            // Call CHPTRD and CUPGTR to compute S and U from AP

            ccopy(NAP, AP, 1, VP, 1 );

            NTEST = 7
            chptrd('L', N, VP, SD, SE, TAU, IINFO );

            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'CHPTRD(L)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 7 ) = ULPINV
                  GO TO 280
               }
            }

            NTEST = 8
            cupgtr('L', N, VP, TAU, U, LDU, WORK, IINFO );
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'CUPGTR(L)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 8 ) = ULPINV
                  GO TO 280
               }
            }

            chpt21(2, 'Lower', N, 1, AP, SD, SE, U, LDU, VP, TAU, WORK, RWORK, RESULT( 7 ) );
            chpt21(3, 'Lower', N, 1, AP, SD, SE, U, LDU, VP, TAU, WORK, RWORK, RESULT( 8 ) );

            // Call CSTEQR to compute D1, D2, and Z, do tests.

            // Compute D1 and Z

            scopy(N, SD, 1, D1, 1 );
            if (N.GT.0) CALL SCOPY( N-1, SE, 1, RWORK, 1 );
            claset('Full', N, N, CZERO, CONE, Z, LDU );

            NTEST = 9
            csteqr('V', N, D1, RWORK, Z, LDU, RWORK( N+1 ), IINFO );
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'CSTEQR(V)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 9 ) = ULPINV
                  GO TO 280
               }
            }

            // Compute D2

            scopy(N, SD, 1, D2, 1 );
            if (N.GT.0) CALL SCOPY( N-1, SE, 1, RWORK, 1 );

            NTEST = 11
            csteqr('N', N, D2, RWORK, WORK, LDU, RWORK( N+1 ), IINFO );
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'CSTEQR(N)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 11 ) = ULPINV
                  GO TO 280
               }
            }

            // Compute D3 (using PWK method)

            scopy(N, SD, 1, D3, 1 );
            if (N.GT.0) CALL SCOPY( N-1, SE, 1, RWORK, 1 );

            NTEST = 12
            ssterf(N, D3, RWORK, IINFO );
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'SSTERF', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 12 ) = ULPINV
                  GO TO 280
               }
            }

            // Do Tests 9 and 10

            cstt21(N, 0, SD, SE, D1, DUMMA, Z, LDU, WORK, RWORK, RESULT( 9 ) );

            // Do Tests 11 and 12

            TEMP1 = ZERO
            TEMP2 = ZERO
            TEMP3 = ZERO
            TEMP4 = ZERO

            for (J = 1; J <= N; J++) { // 150
               TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D2( J ) ) )
               TEMP2 = MAX( TEMP2, ABS( D1( J )-D2( J ) ) )
               TEMP3 = MAX( TEMP3, ABS( D1( J ) ), ABS( D3( J ) ) )
               TEMP4 = MAX( TEMP4, ABS( D1( J )-D3( J ) ) )
            } // 150

            RESULT( 11 ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )
            RESULT( 12 ) = TEMP4 / MAX( UNFL, ULP*MAX( TEMP3, TEMP4 ) )

            // Do Test 13 -- Sturm Sequence Test of Eigenvalues
                          // Go up by factors of two until it succeeds

            NTEST = 13
            TEMP1 = THRESH*( HALF-ULP )

            for (J = 0; J <= LOG2UI; J++) { // 160
               sstech(N, SD, SE, D1, TEMP1, RWORK, IINFO );
               if (IINFO == 0) GO TO 170;
               TEMP1 = TEMP1*TWO
            } // 160

            } // 170
            RESULT( 13 ) = TEMP1

            // For positive definite matrices ( JTYPE.GT.15 ) call CPTEQR
            // and do tests 14, 15, and 16 .

            if ( JTYPE.GT.15 ) {

               // Compute D4 and Z4

               scopy(N, SD, 1, D4, 1 );
               if (N.GT.0) CALL SCOPY( N-1, SE, 1, RWORK, 1 );
               claset('Full', N, N, CZERO, CONE, Z, LDU );

               NTEST = 14
               cpteqr('V', N, D4, RWORK, Z, LDU, RWORK( N+1 ), IINFO );
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'CPTEQR(V)', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( 14 ) = ULPINV
                     GO TO 280
                  }
               }

               // Do Tests 14 and 15

               cstt21(N, 0, SD, SE, D4, DUMMA, Z, LDU, WORK, RWORK, RESULT( 14 ) );

               // Compute D5

               scopy(N, SD, 1, D5, 1 );
               if (N.GT.0) CALL SCOPY( N-1, SE, 1, RWORK, 1 );

               NTEST = 16
               cpteqr('N', N, D5, RWORK, Z, LDU, RWORK( N+1 ), IINFO );
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'CPTEQR(N)', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( 16 ) = ULPINV
                     GO TO 280
                  }
               }

               // Do Test 16

               TEMP1 = ZERO
               TEMP2 = ZERO
               for (J = 1; J <= N; J++) { // 180
                  TEMP1 = MAX( TEMP1, ABS( D4( J ) ), ABS( D5( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( D4( J )-D5( J ) ) )
               } // 180

               RESULT( 16 ) = TEMP2 / MAX( UNFL, HUN*ULP*MAX( TEMP1, TEMP2 ) )
            } else {
               RESULT( 14 ) = ZERO
               RESULT( 15 ) = ZERO
               RESULT( 16 ) = ZERO
            }

            // Call SSTEBZ with different options and do tests 17-18.

               // If S is positive definite and diagonally dominant,
               // ask for all eigenvalues with high relative accuracy.

            VL = ZERO
            VU = ZERO
            IL = 0
            IU = 0
            if ( JTYPE == 21 ) {
               NTEST = 17
               ABSTOL = UNFL + UNFL
               sstebz('A', 'E', N, VL, VU, IL, IU, ABSTOL, SD, SE, M, NSPLIT, WR, IWORK( 1 ), IWORK( N+1 ), RWORK, IWORK( 2*N+1 ), IINFO );
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEBZ(A,rel)', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( 17 ) = ULPINV
                     GO TO 280
                  }
               }

               // Do test 17

               TEMP2 = TWO*( TWO*N-ONE )*ULP*( ONE+EIGHT*HALF**2 ) / ( ONE-HALF )**4

               TEMP1 = ZERO
               for (J = 1; J <= N; J++) { // 190
                  TEMP1 = MAX( TEMP1, ABS( D4( J )-WR( N-J+1 ) ) / ( ABSTOL+ABS( D4( J ) ) ) )
               } // 190

               RESULT( 17 ) = TEMP1 / TEMP2
            } else {
               RESULT( 17 ) = ZERO
            }

            // Now ask for all eigenvalues with high absolute accuracy.

            NTEST = 18
            ABSTOL = UNFL + UNFL
            sstebz('A', 'E', N, VL, VU, IL, IU, ABSTOL, SD, SE, M, NSPLIT, WA1, IWORK( 1 ), IWORK( N+1 ), RWORK, IWORK( 2*N+1 ), IINFO );
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'SSTEBZ(A)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 18 ) = ULPINV
                  GO TO 280
               }
            }

            // Do test 18

            TEMP1 = ZERO
            TEMP2 = ZERO
            for (J = 1; J <= N; J++) { // 200
               TEMP1 = MAX( TEMP1, ABS( D3( J ) ), ABS( WA1( J ) ) )
               TEMP2 = MAX( TEMP2, ABS( D3( J )-WA1( J ) ) )
            } // 200

            RESULT( 18 ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )

            // Choose random values for IL and IU, and ask for the
            // IL-th through IU-th eigenvalues.

            NTEST = 19
            if ( N.LE.1 ) {
               IL = 1
               IU = N
            } else {
               IL = 1 + ( N-1 )*INT( SLARND( 1, ISEED2 ) )
               IU = 1 + ( N-1 )*INT( SLARND( 1, ISEED2 ) )
               if ( IU.LT.IL ) {
                  ITEMP = IU
                  IU = IL
                  IL = ITEMP
               }
            }

            sstebz('I', 'E', N, VL, VU, IL, IU, ABSTOL, SD, SE, M2, NSPLIT, WA2, IWORK( 1 ), IWORK( N+1 ), RWORK, IWORK( 2*N+1 ), IINFO );
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'SSTEBZ(I)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 19 ) = ULPINV
                  GO TO 280
               }
            }

            // Determine the values VL and VU of the IL-th and IU-th
            // eigenvalues and ask for all eigenvalues in this range.

            if ( N.GT.0 ) {
               if ( IL.NE.1 ) {
                  VL = WA1( IL ) - MAX( HALF*( WA1( IL )-WA1( IL-1 ) ), ULP*ANORM, TWO*RTUNFL )
               } else {
                  VL = WA1( 1 ) - MAX( HALF*( WA1( N )-WA1( 1 ) ), ULP*ANORM, TWO*RTUNFL )
               }
               if ( IU.NE.N ) {
                  VU = WA1( IU ) + MAX( HALF*( WA1( IU+1 )-WA1( IU ) ), ULP*ANORM, TWO*RTUNFL )
               } else {
                  VU = WA1( N ) + MAX( HALF*( WA1( N )-WA1( 1 ) ), ULP*ANORM, TWO*RTUNFL )
               }
            } else {
               VL = ZERO
               VU = ONE
            }

            sstebz('V', 'E', N, VL, VU, IL, IU, ABSTOL, SD, SE, M3, NSPLIT, WA3, IWORK( 1 ), IWORK( N+1 ), RWORK, IWORK( 2*N+1 ), IINFO );
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'SSTEBZ(V)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 19 ) = ULPINV
                  GO TO 280
               }
            }

            if ( M3 == 0 .AND. N.NE.0 ) {
               RESULT( 19 ) = ULPINV
               GO TO 280
            }

            // Do test 19

            TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL )
            TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL )
            if ( N.GT.0 ) {
               TEMP3 = MAX( ABS( WA1( N ) ), ABS( WA1( 1 ) ) )
            } else {
               TEMP3 = ZERO
            }

            RESULT( 19 ) = ( TEMP1+TEMP2 ) / MAX( UNFL, TEMP3*ULP )

            // Call CSTEIN to compute eigenvectors corresponding to
            // eigenvalues in WA1.  (First call SSTEBZ again, to make sure
            // it returns these eigenvalues in the correct order.)

            NTEST = 21
            sstebz('A', 'B', N, VL, VU, IL, IU, ABSTOL, SD, SE, M, NSPLIT, WA1, IWORK( 1 ), IWORK( N+1 ), RWORK, IWORK( 2*N+1 ), IINFO );
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'SSTEBZ(A,B)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 20 ) = ULPINV
                  RESULT( 21 ) = ULPINV
                  GO TO 280
               }
            }

            cstein(N, SD, SE, M, WA1, IWORK( 1 ), IWORK( N+1 ), Z, LDU, RWORK, IWORK( 2*N+1 ), IWORK( 3*N+1 ), IINFO );
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'CSTEIN', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 20 ) = ULPINV
                  RESULT( 21 ) = ULPINV
                  GO TO 280
               }
            }

            // Do tests 20 and 21

            cstt21(N, 0, SD, SE, WA1, DUMMA, Z, LDU, WORK, RWORK, RESULT( 20 ) );

            // Call CSTEDC(I) to compute D1 and Z, do tests.

            // Compute D1 and Z

            INDE = 1
            INDRWK = INDE + N
            scopy(N, SD, 1, D1, 1 );
            if (N.GT.0) CALL SCOPY( N-1, SE, 1, RWORK( INDE ), 1 );
            claset('Full', N, N, CZERO, CONE, Z, LDU );

            NTEST = 22
            cstedc('I', N, D1, RWORK( INDE ), Z, LDU, WORK, LWEDC, RWORK( INDRWK ), LRWEDC, IWORK, LIWEDC, IINFO );
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'CSTEDC(I)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 22 ) = ULPINV
                  GO TO 280
               }
            }

            // Do Tests 22 and 23

            cstt21(N, 0, SD, SE, D1, DUMMA, Z, LDU, WORK, RWORK, RESULT( 22 ) );

            // Call CSTEDC(V) to compute D1 and Z, do tests.

            // Compute D1 and Z

            scopy(N, SD, 1, D1, 1 );
            if (N.GT.0) CALL SCOPY( N-1, SE, 1, RWORK( INDE ), 1 );
            claset('Full', N, N, CZERO, CONE, Z, LDU );

            NTEST = 24
            cstedc('V', N, D1, RWORK( INDE ), Z, LDU, WORK, LWEDC, RWORK( INDRWK ), LRWEDC, IWORK, LIWEDC, IINFO );
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'CSTEDC(V)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 24 ) = ULPINV
                  GO TO 280
               }
            }

            // Do Tests 24 and 25

            cstt21(N, 0, SD, SE, D1, DUMMA, Z, LDU, WORK, RWORK, RESULT( 24 ) );

            // Call CSTEDC(N) to compute D2, do tests.

            // Compute D2

            scopy(N, SD, 1, D2, 1 );
            if (N.GT.0) CALL SCOPY( N-1, SE, 1, RWORK( INDE ), 1 );
            claset('Full', N, N, CZERO, CONE, Z, LDU );

            NTEST = 26
            cstedc('N', N, D2, RWORK( INDE ), Z, LDU, WORK, LWEDC, RWORK( INDRWK ), LRWEDC, IWORK, LIWEDC, IINFO );
            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'CSTEDC(N)', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               if ( IINFO.LT.0 ) {
                  RETURN
               } else {
                  RESULT( 26 ) = ULPINV
                  GO TO 280
               }
            }

            // Do Test 26

            TEMP1 = ZERO
            TEMP2 = ZERO

            for (J = 1; J <= N; J++) { // 210
               TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D2( J ) ) )
               TEMP2 = MAX( TEMP2, ABS( D1( J )-D2( J ) ) )
            } // 210

            RESULT( 26 ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )

            // Only test CSTEMR if IEEE compliant

            if ( ILAENV( 10, 'CSTEMR', 'VA', 1, 0, 0, 0 ) == 1 .AND. ILAENV( 11, 'CSTEMR', 'VA', 1, 0, 0, 0 ) == 1 ) {

            // Call CSTEMR, do test 27 (relative eigenvalue accuracy)

               // If S is positive definite and diagonally dominant,
               // ask for all eigenvalues with high relative accuracy.

               VL = ZERO
               VU = ZERO
               IL = 0
               IU = 0
               if ( JTYPE == 21 .AND. CREL ) {
                  NTEST = 27
                  ABSTOL = UNFL + UNFL
                  cstemr('V', 'A', N, SD, SE, VL, VU, IL, IU, M, WR, Z, LDU, N, IWORK( 1 ), TRYRAC, RWORK, LRWORK, IWORK( 2*N+1 ), LWORK-2*N, IINFO );
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'CSTEMR(V,A,rel)', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( 27 ) = ULPINV
                        GO TO 270
                     }
                  }

               // Do test 27

                  TEMP2 = TWO*( TWO*N-ONE )*ULP*( ONE+EIGHT*HALF**2 ) / ( ONE-HALF )**4

                  TEMP1 = ZERO
                  for (J = 1; J <= N; J++) { // 220
                     TEMP1 = MAX( TEMP1, ABS( D4( J )-WR( N-J+1 ) ) / ( ABSTOL+ABS( D4( J ) ) ) )
                  } // 220

                  RESULT( 27 ) = TEMP1 / TEMP2

                  IL = 1 + ( N-1 )*INT( SLARND( 1, ISEED2 ) )
                  IU = 1 + ( N-1 )*INT( SLARND( 1, ISEED2 ) )
                  if ( IU.LT.IL ) {
                     ITEMP = IU
                     IU = IL
                     IL = ITEMP
                  }

                  if ( CRANGE ) {
                     NTEST = 28
                     ABSTOL = UNFL + UNFL
                     cstemr('V', 'I', N, SD, SE, VL, VU, IL, IU, M, WR, Z, LDU, N, IWORK( 1 ), TRYRAC, RWORK, LRWORK, IWORK( 2*N+1 ), LWORK-2*N, IINFO );

                     if ( IINFO.NE.0 ) {
                        WRITE( NOUNIT, FMT = 9999 )'CSTEMR(V,I,rel)', IINFO, N, JTYPE, IOLDSD
                        INFO = ABS( IINFO )
                        if ( IINFO.LT.0 ) {
                           RETURN
                        } else {
                           RESULT( 28 ) = ULPINV
                           GO TO 270
                        }
                     }


                  // Do test 28

                     TEMP2 = TWO*( TWO*N-ONE )*ULP* ( ONE+EIGHT*HALF**2 ) / ( ONE-HALF )**4

                     TEMP1 = ZERO
                     for (J = IL; J <= IU; J++) { // 230
                        TEMP1 = MAX( TEMP1, ABS( WR( J-IL+1 )-D4( N-J+ 1 ) ) / ( ABSTOL+ABS( WR( J-IL+1 ) ) ) )
                     } // 230

                     RESULT( 28 ) = TEMP1 / TEMP2
                  } else {
                     RESULT( 28 ) = ZERO
                  }
               } else {
                  RESULT( 27 ) = ZERO
                  RESULT( 28 ) = ZERO
               }

            // Call CSTEMR(V,I) to compute D1 and Z, do tests.

            // Compute D1 and Z

               scopy(N, SD, 1, D5, 1 );
               if (N.GT.0) CALL SCOPY( N-1, SE, 1, RWORK, 1 );
               claset('Full', N, N, CZERO, CONE, Z, LDU );

               if ( CRANGE ) {
                  NTEST = 29
                  IL = 1 + ( N-1 )*INT( SLARND( 1, ISEED2 ) )
                  IU = 1 + ( N-1 )*INT( SLARND( 1, ISEED2 ) )
                  if ( IU.LT.IL ) {
                     ITEMP = IU
                     IU = IL
                     IL = ITEMP
                  }
                  cstemr('V', 'I', N, D5, RWORK, VL, VU, IL, IU, M, D1, Z, LDU, N, IWORK( 1 ), TRYRAC, RWORK( N+1 ), LRWORK-N, IWORK( 2*N+1 ), LIWORK-2*N, IINFO );
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'CSTEMR(V,I)', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( 29 ) = ULPINV
                        GO TO 280
                     }
                  }

            // Do Tests 29 and 30


            // Call CSTEMR to compute D2, do tests.

            // Compute D2

                  scopy(N, SD, 1, D5, 1 );
                  if (N.GT.0) CALL SCOPY( N-1, SE, 1, RWORK, 1 );

                  NTEST = 31
                  cstemr('N', 'I', N, D5, RWORK, VL, VU, IL, IU, M, D2, Z, LDU, N, IWORK( 1 ), TRYRAC, RWORK( N+1 ), LRWORK-N, IWORK( 2*N+1 ), LIWORK-2*N, IINFO );
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'CSTEMR(N,I)', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( 31 ) = ULPINV
                        GO TO 280
                     }
                  }

            // Do Test 31

                  TEMP1 = ZERO
                  TEMP2 = ZERO

                  for (J = 1; J <= IU - IL + 1; J++) { // 240
                     TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D2( J ) ) )
                     TEMP2 = MAX( TEMP2, ABS( D1( J )-D2( J ) ) )
                  } // 240

                  RESULT( 31 ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )


            // Call CSTEMR(V,V) to compute D1 and Z, do tests.

            // Compute D1 and Z

                  scopy(N, SD, 1, D5, 1 );
                  if (N.GT.0) CALL SCOPY( N-1, SE, 1, RWORK, 1 );
                  claset('Full', N, N, CZERO, CONE, Z, LDU );

                  NTEST = 32

                  if ( N.GT.0 ) {
                     if ( IL.NE.1 ) {
                        VL = D2( IL ) - MAX( HALF* ( D2( IL )-D2( IL-1 ) ), ULP*ANORM, TWO*RTUNFL )
                     } else {
                        VL = D2( 1 ) - MAX( HALF*( D2( N )-D2( 1 ) ), ULP*ANORM, TWO*RTUNFL )
                     }
                     if ( IU.NE.N ) {
                        VU = D2( IU ) + MAX( HALF* ( D2( IU+1 )-D2( IU ) ), ULP*ANORM, TWO*RTUNFL )
                     } else {
                        VU = D2( N ) + MAX( HALF*( D2( N )-D2( 1 ) ), ULP*ANORM, TWO*RTUNFL )
                     }
                  } else {
                     VL = ZERO
                     VU = ONE
                  }

                  cstemr('V', 'V', N, D5, RWORK, VL, VU, IL, IU, M, D1, Z, LDU, N, IWORK( 1 ), TRYRAC, RWORK( N+1 ), LRWORK-N, IWORK( 2*N+1 ), LIWORK-2*N, IINFO );
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'CSTEMR(V,V)', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( 32 ) = ULPINV
                        GO TO 280
                     }
                  }

            // Do Tests 32 and 33

                  cstt22(N, M, 0, SD, SE, D1, DUMMA, Z, LDU, WORK, M, RWORK, RESULT( 32 ) );

            // Call CSTEMR to compute D2, do tests.

            // Compute D2

                  scopy(N, SD, 1, D5, 1 );
                  if (N.GT.0) CALL SCOPY( N-1, SE, 1, RWORK, 1 );

                  NTEST = 34
                  cstemr('N', 'V', N, D5, RWORK, VL, VU, IL, IU, M, D2, Z, LDU, N, IWORK( 1 ), TRYRAC, RWORK( N+1 ), LRWORK-N, IWORK( 2*N+1 ), LIWORK-2*N, IINFO );
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'CSTEMR(N,V)', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( 34 ) = ULPINV
                        GO TO 280
                     }
                  }

            // Do Test 34

                  TEMP1 = ZERO
                  TEMP2 = ZERO

                  for (J = 1; J <= IU - IL + 1; J++) { // 250
                     TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D2( J ) ) )
                     TEMP2 = MAX( TEMP2, ABS( D1( J )-D2( J ) ) )
                  } // 250

                  RESULT( 34 ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )
               } else {
                  RESULT( 29 ) = ZERO
                  RESULT( 30 ) = ZERO
                  RESULT( 31 ) = ZERO
                  RESULT( 32 ) = ZERO
                  RESULT( 33 ) = ZERO
                  RESULT( 34 ) = ZERO
               }


            // Call CSTEMR(V,A) to compute D1 and Z, do tests.

            // Compute D1 and Z

               scopy(N, SD, 1, D5, 1 );
               if (N.GT.0) CALL SCOPY( N-1, SE, 1, RWORK, 1 );

               NTEST = 35

               cstemr('V', 'A', N, D5, RWORK, VL, VU, IL, IU, M, D1, Z, LDU, N, IWORK( 1 ), TRYRAC, RWORK( N+1 ), LRWORK-N, IWORK( 2*N+1 ), LIWORK-2*N, IINFO );
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'CSTEMR(V,A)', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( 35 ) = ULPINV
                     GO TO 280
                  }
               }

            // Do Tests 35 and 36

               cstt22(N, M, 0, SD, SE, D1, DUMMA, Z, LDU, WORK, M, RWORK, RESULT( 35 ) );

            // Call CSTEMR to compute D2, do tests.

            // Compute D2

               scopy(N, SD, 1, D5, 1 );
               if (N.GT.0) CALL SCOPY( N-1, SE, 1, RWORK, 1 );

               NTEST = 37
               cstemr('N', 'A', N, D5, RWORK, VL, VU, IL, IU, M, D2, Z, LDU, N, IWORK( 1 ), TRYRAC, RWORK( N+1 ), LRWORK-N, IWORK( 2*N+1 ), LIWORK-2*N, IINFO );
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'CSTEMR(N,A)', IINFO, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  if ( IINFO.LT.0 ) {
                     RETURN
                  } else {
                     RESULT( 37 ) = ULPINV
                     GO TO 280
                  }
               }

            // Do Test 34

               TEMP1 = ZERO
               TEMP2 = ZERO

               for (J = 1; J <= N; J++) { // 260
                  TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D2( J ) ) )
                  TEMP2 = MAX( TEMP2, ABS( D1( J )-D2( J ) ) )
               } // 260

               RESULT( 37 ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )
            }
            } // 270
            } // 280
            NTESTT = NTESTT + NTEST

            // End of Loop -- Check for RESULT(j) > THRESH


            // Print out tests which fail.

            for (JR = 1; JR <= NTEST; JR++) { // 290
               if ( RESULT( JR ).GE.THRESH ) {

                  // If this is the first test to fail,
                  // print a header to the data file.

                  if ( NERRS == 0 ) {
                     WRITE( NOUNIT, FMT = 9998 )'CST'
                     WRITE( NOUNIT, FMT = 9997 )
                     WRITE( NOUNIT, FMT = 9996 )
                     WRITE( NOUNIT, FMT = 9995 )'Hermitian'
                     WRITE( NOUNIT, FMT = 9994 )

                     // Tests performed

                     WRITE( NOUNIT, FMT = 9987 )
                  }
                  NERRS = NERRS + 1
                  if ( RESULT( JR ).LT.10000.0E0 ) {
                     WRITE( NOUNIT, FMT = 9989 )N, JTYPE, IOLDSD, JR, RESULT( JR )
                  } else {
                     WRITE( NOUNIT, FMT = 9988 )N, JTYPE, IOLDSD, JR, RESULT( JR )
                  }
               }
            } // 290
         } // 300
      } // 310

      // Summary

      slasum('CST', NOUNIT, NERRS, NTESTT );
      RETURN

 9999 FORMAT( ' CCHKST: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )

 9998 FORMAT( / 1X, A3, ' -- Complex Hermitian eigenvalue problem' )
 9997 FORMAT( ' Matrix types (see CCHKST for details): ' )

 9996 FORMAT( / ' Special Matrices:', / '  1=Zero matrix.                        ', '  5=Diagonal: clustered entries.', / '  2=Identity matrix.                    ', '  6=Diagonal: large, evenly spaced.', / '  3=Diagonal: evenly spaced entries.    ', '  7=Diagonal: small, evenly spaced.', / '  4=Diagonal: geometr. spaced entries.' )
 9995 FORMAT( ' Dense ', A, ' Matrices:', / '  8=Evenly spaced eigenvals.            ', ' 12=Small, evenly spaced eigenvals.', / '  9=Geometrically spaced eigenvals.     ', ' 13=Matrix with random O(1) entries.', / ' 10=Clustered eigenvalues.              ', ' 14=Matrix with large random entries.', / ' 11=Large, evenly spaced eigenvals.     ', ' 15=Matrix with small random entries.' )
 9994 FORMAT( ' 16=Positive definite, evenly spaced eigenvalues', / ' 17=Positive definite, geometrically spaced eigenvlaues', / ' 18=Positive definite, clustered eigenvalues', / ' 19=Positive definite, small evenly spaced eigenvalues', / ' 20=Positive definite, large evenly spaced eigenvalues', / ' 21=Diagonally dominant tridiagonal, geometrically', ' spaced eigenvalues' )

 9989 FORMAT( ' Matrix order=', I5, ', type=', I2, ', seed=', 4( I4, ',' ), ' result ', I3, ' is', 0P, F8.2 )
 9988 FORMAT( ' Matrix order=', I5, ', type=', I2, ', seed=', 4( I4, ',' ), ' result ', I3, ' is', 1P, E10.3 )

 9987 FORMAT( / 'Test performed:  see CCHKST for details.', / )
      // End of CCHKST

      }
