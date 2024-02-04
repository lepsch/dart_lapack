      void dchkst2stg(NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH, NOUNIT, A, LDA, AP, SD, SE, D1, D2, D3, D4, D5, WA1, WA2, WA3, WR, U, LDU, V, VP, TAU, Z, WORK, LWORK, IWORK, LIWORK, RESULT, INFO ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDU, LIWORK, LWORK, NOUNIT, NSIZES, NTYPES;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                ISEED( 4 ), IWORK( * ), NN( * );
      double             A( LDA, * ), AP( * ), D1( * ), D2( * ), D3( * ), D4( * ), D5( * ), RESULT( * ), SD( * ), SE( * ), TAU( * ), U( LDU, * ), V( LDU, * ), VP( * ), WA1( * ), WA2( * ), WA3( * ), WORK( * ), WR( * ), Z( LDU, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TWO, EIGHT, TEN, HUN;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0, EIGHT = 8.0, TEN = 10.0, HUN = 100.0 ;
      double             HALF;
      const              HALF = ONE / TWO ;
      int                MAXTYP;
      const              MAXTYP = 21 ;
      bool               SRANGE;
      const              SRANGE = false ;
      bool               SREL;
      const              SREL = false ;
      // ..
      // .. Local Scalars ..
      bool               BADNN, TRYRAC;
      int                I, IINFO, IL, IMODE, ITEMP, ITYPE, IU, J, JC, JR, JSIZE, JTYPE, LGN, LIWEDC, LOG2UI, LWEDC, M, M2, M3, MTYPES, N, NAP, NBLOCK, NERRS, NMATS, NMAX, NSPLIT, NTEST, NTESTT, LH, LW;
      double             ABSTOL, ANINV, ANORM, COND, OVFL, RTOVFL, RTUNFL, TEMP1, TEMP2, TEMP3, TEMP4, ULP, ULPINV, UNFL, VL, VU;
      // ..
      // .. Local Arrays ..
      int                IDUMMA( 1 ), IOLDSD( 4 ), ISEED2( 4 ), KMAGN( MAXTYP ), KMODE( MAXTYP ), KTYPE( MAXTYP );
      double             DUMMA( 1 );
      // ..
      // .. External Functions ..
      //- int                ILAENV;
      //- double             DLAMCH, DLARND, DSXT1;
      // EXTERNAL ILAENV, DLAMCH, DLARND, DSXT1
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DLACPY, DLASET, DLASUM, DLATMR, DLATMS, DOPGTR, DORGTR, DPTEQR, DSPT21, DSPTRD, DSTEBZ, DSTECH, DSTEDC, DSTEMR, DSTEIN, DSTEQR, DSTERF, DSTT21, DSTT22, DSYT21, DSYTRD, XERBLA, DSYTRD_2STAGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, INT, LOG, MAX, MIN, SQRT
      // ..
      // .. Data statements ..
      const KTYPE = [ 1, 2, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 8, 8, 8, 9, 9, 9, 9, 9, 10 ];
      const KMAGN = [ 1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3, 1, 1, 1, 2, 3, 1 ];
      const KMODE = [ 0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0, 4, 3, 1, 4, 4, 3 ];
      // ..
      // .. Executable Statements ..

      // Keep ftnchek happy
      IDUMMA[1] = 1;

      // Check for errors

      NTESTT = 0;
      INFO = 0;

      // Important constants

      BADNN = false;
      TRYRAC = true;
      NMAX = 1;
      for (J = 1; J <= NSIZES; J++) { // 10
         NMAX = max( NMAX, NN( J ) );
         if( NN( J ) < 0 ) BADNN = true;
      } // 10

      NBLOCK = ILAENV( 1, 'DSYTRD', 'L', NMAX, -1, -1, -1 );
      NBLOCK = min( NMAX, max( 1, NBLOCK ) );

      // Check for errors

      if ( NSIZES < 0 ) {
         INFO = -1;
      } else if ( BADNN ) {
         INFO = -2;
      } else if ( NTYPES < 0 ) {
         INFO = -3;
      } else if ( LDA < NMAX ) {
         INFO = -9;
      } else if ( LDU < NMAX ) {
         INFO = -23;
      } else if ( 2*max( 2, NMAX )**2 > LWORK ) {
         INFO = -29;
      }

      if ( INFO != 0 ) {
         xerbla('DCHKST2STG', -INFO );
         return;
      }

      // Quick return if possible

      if (NSIZES == 0 || NTYPES == 0) return;

      // More Important constants

      UNFL = DLAMCH( 'Safe minimum' );
      OVFL = ONE / UNFL;
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' );
      ULPINV = ONE / ULP;
      LOG2UI = INT( LOG( ULPINV ) / LOG( TWO ) );
      RTUNFL = sqrt( UNFL );
      RTOVFL = sqrt( OVFL );

      // Loop over sizes, types

      for (I = 1; I <= 4; I++) { // 20
         ISEED2[I] = ISEED( I );
      } // 20
      NERRS = 0;
      NMATS = 0;

      for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) { // 310
         N = NN( JSIZE );
         if ( N > 0 ) {
            LGN = INT( LOG( N.toDouble() ) / LOG( TWO ) );
            if (2**LGN < N) LGN = LGN + 1;
            IF( 2**LGN < N ) LGN = LGN + 1;
            LWEDC = 1 + 4*N + 2*N*LGN + 4*N**2;
            LIWEDC = 6 + 6*N + 5*N*LGN;
         } else {
            LWEDC = 8;
            LIWEDC = 12;
         }
         NAP = ( N*( N+1 ) ) / 2;
         ANINV = ONE / (max( 1, N )).toDouble();

         if ( NSIZES != 1 ) {
            MTYPES = min( MAXTYP, NTYPES );
         } else {
            MTYPES = min( MAXTYP+1, NTYPES );
         }

         for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) { // 300
            if( !DOTYPE( JTYPE ) ) GO TO 300;
            NMATS = NMATS + 1;
            NTEST = 0;

            for (J = 1; J <= 4; J++) { // 30
               IOLDSD[J] = ISEED( J );
            } // 30

            // Compute "A"

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

            if (MTYPES > MAXTYP) GO TO 100;

            ITYPE = KTYPE( JTYPE );
            IMODE = KMODE( JTYPE );

            // Compute norm

            GO TO ( 40, 50, 60 )KMAGN( JTYPE );

            } // 40
            ANORM = ONE;
            GO TO 70;

            } // 50
            ANORM = ( RTOVFL*ULP )*ANINV;
            GO TO 70;

            } // 60
            ANORM = RTUNFL*N*ULPINV;
            GO TO 70;

            } // 70

            dlaset('Full', LDA, N, ZERO, ZERO, A, LDA );
            IINFO = 0;
            if ( JTYPE <= 15 ) {
               COND = ULPINV;
            } else {
               COND = ULPINV*ANINV / TEN;
            }

            // Special Matrices -- Identity & Jordan block

               // Zero

            if ( ITYPE == 1 ) {
               IINFO = 0;

            } else if ( ITYPE == 2 ) {

               // Identity

               for (JC = 1; JC <= N; JC++) { // 80
                  A[JC, JC] = ANORM;
               } // 80

            } else if ( ITYPE == 4 ) {

               // Diagonal Matrix, [Eigen]values Specified

               dlatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, 0, 0, 'N', A, LDA, WORK( N+1 ), IINFO );


            } else if ( ITYPE == 5 ) {

               // Symmetric, eigenvalues specified

               dlatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, N, N, 'N', A, LDA, WORK( N+1 ), IINFO );

            } else if ( ITYPE == 7 ) {

               // Diagonal, random eigenvalues

               dlatmr(N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, 0, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else if ( ITYPE == 8 ) {

               // Symmetric, random eigenvalues

               dlatmr(N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else if ( ITYPE == 9 ) {

               // Positive definite, eigenvalues specified.

               dlatms(N, N, 'S', ISEED, 'P', WORK, IMODE, COND, ANORM, N, N, 'N', A, LDA, WORK( N+1 ), IINFO );

            } else if ( ITYPE == 10 ) {

               // Positive definite tridiagonal, eigenvalues specified.

               dlatms(N, N, 'S', ISEED, 'P', WORK, IMODE, COND, ANORM, 1, 1, 'N', A, LDA, WORK( N+1 ), IINFO );
               for (I = 2; I <= N; I++) { // 90
                  TEMP1 = ( A( I-1, I ) ).abs() / sqrt( ABS( A( I-1, I-1 )*A( I, I ) ) );
                  if ( TEMP1 > HALF ) {
                     A[I-1, I] = HALF*sqrt( ABS( A( I-1, I-1 )*A( I, I ) ) );
                     A[I, I-1] = A( I-1, I );
                  }
               } // 90

            } else {

               IINFO = 1;
            }

            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'Generator', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               return;
            }

            } // 100

            // Call DSYTRD and DORGTR to compute S and U from
            // upper triangle.

            dlacpy('U', N, N, A, LDA, V, LDU );

            NTEST = 1;
            dsytrd('U', N, V, LDU, SD, SE, TAU, WORK, LWORK, IINFO );

            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'DSYTRD(U)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               if ( IINFO < 0 ) {
                  return;
               } else {
                  RESULT[1] = ULPINV;
                  GO TO 280;
               }
            }

            dlacpy('U', N, N, V, LDU, U, LDU );

            NTEST = 2;
            dorgtr('U', N, U, LDU, TAU, WORK, LWORK, IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'DORGTR(U)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               if ( IINFO < 0 ) {
                  return;
               } else {
                  RESULT[2] = ULPINV;
                  GO TO 280;
               }
            }

            // Do tests 1 and 2

            dsyt21(2, 'Upper', N, 1, A, LDA, SD, SE, U, LDU, V, LDU, TAU, WORK, RESULT( 1 ) );
            dsyt21(3, 'Upper', N, 1, A, LDA, SD, SE, U, LDU, V, LDU, TAU, WORK, RESULT( 2 ) );

            // Compute D1 the eigenvalues resulting from the tridiagonal
            // form using the standard 1-stage algorithm and use it as a
            // reference to compare with the 2-stage technique

            // Compute D1 from the 1-stage and used as reference for the
            // 2-stage

            dcopy(N, SD, 1, D1, 1 );
            if (N > 0) dcopy( N-1, SE, 1, WORK, 1 );

            dsteqr('N', N, D1, WORK, WORK( N+1 ), LDU, WORK( N+1 ), IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'DSTEQR(N)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               if ( IINFO < 0 ) {
                  return;
               } else {
                  RESULT[3] = ULPINV;
                  GO TO 280;
               }
            }

            // 2-STAGE TRD Upper case is used to compute D2.
            // Note to set SD and SE to zero to be sure not reusing
            // the one from above. Compare it with D1 computed
            // using the 1-stage.

            dlaset('Full', N, 1, ZERO, ZERO, SD, N );
            dlaset('Full', N, 1, ZERO, ZERO, SE, N );
            dlacpy("U", N, N, A, LDA, V, LDU );
            LH = max(1, 4*N);
            LW = LWORK - LH;
            dsytrd_2stage('N', "U", N, V, LDU, SD, SE, TAU,  WORK, LH, WORK( LH+1 ), LW, IINFO );

            // Compute D2 from the 2-stage Upper case

            dcopy(N, SD, 1, D2, 1 );
            if (N > 0) dcopy( N-1, SE, 1, WORK, 1 );

            dsteqr('N', N, D2, WORK, WORK( N+1 ), LDU, WORK( N+1 ), IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'DSTEQR(N)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               if ( IINFO < 0 ) {
                  return;
               } else {
                  RESULT[3] = ULPINV;
                  GO TO 280;
               }
            }

            // 2-STAGE TRD Lower case is used to compute D3.
            // Note to set SD and SE to zero to be sure not reusing
            // the one from above. Compare it with D1 computed
            // using the 1-stage.

            dlaset('Full', N, 1, ZERO, ZERO, SD, N );
            dlaset('Full', N, 1, ZERO, ZERO, SE, N );
            dlacpy("L", N, N, A, LDA, V, LDU );
            dsytrd_2stage('N', "L", N, V, LDU, SD, SE, TAU,  WORK, LH, WORK( LH+1 ), LW, IINFO );

            // Compute D3 from the 2-stage Upper case

            dcopy(N, SD, 1, D3, 1 );
            if (N > 0) dcopy( N-1, SE, 1, WORK, 1 );

            dsteqr('N', N, D3, WORK, WORK( N+1 ), LDU, WORK( N+1 ), IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'DSTEQR(N)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               if ( IINFO < 0 ) {
                  return;
               } else {
                  RESULT[4] = ULPINV;
                  GO TO 280;
               }
            }

            // Do Tests 3 and 4 which are similar to 11 and 12 but with the
            // D1 computed using the standard 1-stage reduction as reference

            NTEST = 4;
            TEMP1 = ZERO;
            TEMP2 = ZERO;
            TEMP3 = ZERO;
            TEMP4 = ZERO;

            for (J = 1; J <= N; J++) { // 151
               TEMP1 = max( TEMP1, ( D1( J ) ).abs(), ( D2( J ) ) ).abs();
               TEMP2 = max( TEMP2, ABS( D1( J )-D2( J ) ) );
               TEMP3 = max( TEMP3, ( D1( J ) ).abs(), ( D3( J ) ) ).abs();
               TEMP4 = max( TEMP4, ABS( D1( J )-D3( J ) ) );
            } // 151

            RESULT[3] = TEMP2 / max( UNFL, ULP*max( TEMP1, TEMP2 ) );
            RESULT[4] = TEMP4 / max( UNFL, ULP*max( TEMP3, TEMP4 ) );

            // Store the upper triangle of A in AP

            I = 0;
            for (JC = 1; JC <= N; JC++) { // 120
               for (JR = 1; JR <= JC; JR++) { // 110
                  I = I + 1;
                  AP[I] = A( JR, JC );
               } // 110
            } // 120

            // Call DSPTRD and DOPGTR to compute S and U from AP

            dcopy(NAP, AP, 1, VP, 1 );

            NTEST = 5;
            dsptrd('U', N, VP, SD, SE, TAU, IINFO );

            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'DSPTRD(U)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               if ( IINFO < 0 ) {
                  return;
               } else {
                  RESULT[5] = ULPINV;
                  GO TO 280;
               }
            }

            NTEST = 6;
            dopgtr('U', N, VP, TAU, U, LDU, WORK, IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'DOPGTR(U)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               if ( IINFO < 0 ) {
                  return;
               } else {
                  RESULT[6] = ULPINV;
                  GO TO 280;
               }
            }

            // Do tests 5 and 6

            dspt21(2, 'Upper', N, 1, AP, SD, SE, U, LDU, VP, TAU, WORK, RESULT( 5 ) );
            dspt21(3, 'Upper', N, 1, AP, SD, SE, U, LDU, VP, TAU, WORK, RESULT( 6 ) );

            // Store the lower triangle of A in AP

            I = 0;
            for (JC = 1; JC <= N; JC++) { // 140
               for (JR = JC; JR <= N; JR++) { // 130
                  I = I + 1;
                  AP[I] = A( JR, JC );
               } // 130
            } // 140

            // Call DSPTRD and DOPGTR to compute S and U from AP

            dcopy(NAP, AP, 1, VP, 1 );

            NTEST = 7;
            dsptrd('L', N, VP, SD, SE, TAU, IINFO );

            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'DSPTRD(L)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               if ( IINFO < 0 ) {
                  return;
               } else {
                  RESULT[7] = ULPINV;
                  GO TO 280;
               }
            }

            NTEST = 8;
            dopgtr('L', N, VP, TAU, U, LDU, WORK, IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'DOPGTR(L)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               if ( IINFO < 0 ) {
                  return;
               } else {
                  RESULT[8] = ULPINV;
                  GO TO 280;
               }
            }

            dspt21(2, 'Lower', N, 1, AP, SD, SE, U, LDU, VP, TAU, WORK, RESULT( 7 ) );
            dspt21(3, 'Lower', N, 1, AP, SD, SE, U, LDU, VP, TAU, WORK, RESULT( 8 ) );

            // Call DSTEQR to compute D1, D2, and Z, do tests.

            // Compute D1 and Z

            dcopy(N, SD, 1, D1, 1 );
            if (N > 0) dcopy( N-1, SE, 1, WORK, 1 );
            dlaset('Full', N, N, ZERO, ONE, Z, LDU );

            NTEST = 9;
            dsteqr('V', N, D1, WORK, Z, LDU, WORK( N+1 ), IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'DSTEQR(V)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               if ( IINFO < 0 ) {
                  return;
               } else {
                  RESULT[9] = ULPINV;
                  GO TO 280;
               }
            }

            // Compute D2

            dcopy(N, SD, 1, D2, 1 );
            if (N > 0) dcopy( N-1, SE, 1, WORK, 1 );

            NTEST = 11;
            dsteqr('N', N, D2, WORK, WORK( N+1 ), LDU, WORK( N+1 ), IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'DSTEQR(N)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               if ( IINFO < 0 ) {
                  return;
               } else {
                  RESULT[11] = ULPINV;
                  GO TO 280;
               }
            }

            // Compute D3 (using PWK method)

            dcopy(N, SD, 1, D3, 1 );
            if (N > 0) dcopy( N-1, SE, 1, WORK, 1 );

            NTEST = 12;
            dsterf(N, D3, WORK, IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'DSTERF', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               if ( IINFO < 0 ) {
                  return;
               } else {
                  RESULT[12] = ULPINV;
                  GO TO 280;
               }
            }

            // Do Tests 9 and 10

            dstt21(N, 0, SD, SE, D1, DUMMA, Z, LDU, WORK, RESULT( 9 ) );

            // Do Tests 11 and 12

            TEMP1 = ZERO;
            TEMP2 = ZERO;
            TEMP3 = ZERO;
            TEMP4 = ZERO;

            for (J = 1; J <= N; J++) { // 150
               TEMP1 = max( TEMP1, ( D1( J ) ).abs(), ( D2( J ) ) ).abs();
               TEMP2 = max( TEMP2, ABS( D1( J )-D2( J ) ) );
               TEMP3 = max( TEMP3, ( D1( J ) ).abs(), ( D3( J ) ) ).abs();
               TEMP4 = max( TEMP4, ABS( D1( J )-D3( J ) ) );
            } // 150

            RESULT[11] = TEMP2 / max( UNFL, ULP*max( TEMP1, TEMP2 ) );
            RESULT[12] = TEMP4 / max( UNFL, ULP*max( TEMP3, TEMP4 ) );

            // Do Test 13 -- Sturm Sequence Test of Eigenvalues
                          // Go up by factors of two until it succeeds

            NTEST = 13;
            TEMP1 = THRESH*( HALF-ULP );

            for (J = 0; J <= LOG2UI; J++) { // 160
               dstech(N, SD, SE, D1, TEMP1, WORK, IINFO );
               if (IINFO == 0) GO TO 170;
               TEMP1 = TEMP1*TWO;
            } // 160

            } // 170
            RESULT[13] = TEMP1;

            // For positive definite matrices ( JTYPE > 15 ) call DPTEQR
            // and do tests 14, 15, and 16 .

            if ( JTYPE > 15 ) {

               // Compute D4 and Z4

               dcopy(N, SD, 1, D4, 1 );
               if (N > 0) dcopy( N-1, SE, 1, WORK, 1 );
               dlaset('Full', N, N, ZERO, ONE, Z, LDU );

               NTEST = 14;
               dpteqr('V', N, D4, WORK, Z, LDU, WORK( N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'DPTEQR(V)', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[14] = ULPINV;
                     GO TO 280;
                  }
               }

               // Do Tests 14 and 15

               dstt21(N, 0, SD, SE, D4, DUMMA, Z, LDU, WORK, RESULT( 14 ) );

               // Compute D5

               dcopy(N, SD, 1, D5, 1 );
               if (N > 0) dcopy( N-1, SE, 1, WORK, 1 );

               NTEST = 16;
               dpteqr('N', N, D5, WORK, Z, LDU, WORK( N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'DPTEQR(N)', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[16] = ULPINV;
                     GO TO 280;
                  }
               }

               // Do Test 16

               TEMP1 = ZERO;
               TEMP2 = ZERO;
               for (J = 1; J <= N; J++) { // 180
                  TEMP1 = max( TEMP1, ( D4( J ) ).abs(), ( D5( J ) ) ).abs();
                  TEMP2 = max( TEMP2, ABS( D4( J )-D5( J ) ) );
               } // 180

               RESULT[16] = TEMP2 / max( UNFL, HUN*ULP*max( TEMP1, TEMP2 ) );
            } else {
               RESULT[14] = ZERO;
               RESULT[15] = ZERO;
               RESULT[16] = ZERO;
            }

            // Call DSTEBZ with different options and do tests 17-18.

               // If S is positive definite and diagonally dominant,
               // ask for all eigenvalues with high relative accuracy.

            VL = ZERO;
            VU = ZERO;
            IL = 0;
            IU = 0;
            if ( JTYPE == 21 ) {
               NTEST = 17;
               ABSTOL = UNFL + UNFL;
               dstebz('A', 'E', N, VL, VU, IL, IU, ABSTOL, SD, SE, M, NSPLIT, WR, IWORK( 1 ), IWORK( N+1 ), WORK, IWORK( 2*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'DSTEBZ(A,rel)', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[17] = ULPINV;
                     GO TO 280;
                  }
               }

               // Do test 17

               TEMP2 = TWO*( TWO*N-ONE )*ULP*( ONE+EIGHT*HALF**2 ) / ( ONE-HALF )**4;

               TEMP1 = ZERO;
               for (J = 1; J <= N; J++) { // 190
                  TEMP1 = max( TEMP1, ABS( D4( J )-WR( N-J+1 ) ) / ( ABSTOL+( D4( J ) ) ) ).abs();
               } // 190

               RESULT[17] = TEMP1 / TEMP2;
            } else {
               RESULT[17] = ZERO;
            }

            // Now ask for all eigenvalues with high absolute accuracy.

            NTEST = 18;
            ABSTOL = UNFL + UNFL;
            dstebz('A', 'E', N, VL, VU, IL, IU, ABSTOL, SD, SE, M, NSPLIT, WA1, IWORK( 1 ), IWORK( N+1 ), WORK, IWORK( 2*N+1 ), IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'DSTEBZ(A)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               if ( IINFO < 0 ) {
                  return;
               } else {
                  RESULT[18] = ULPINV;
                  GO TO 280;
               }
            }

            // Do test 18

            TEMP1 = ZERO;
            TEMP2 = ZERO;
            for (J = 1; J <= N; J++) { // 200
               TEMP1 = max( TEMP1, ( D3( J ) ).abs(), ( WA1( J ) ) ).abs();
               TEMP2 = max( TEMP2, ABS( D3( J )-WA1( J ) ) );
            } // 200

            RESULT[18] = TEMP2 / max( UNFL, ULP*max( TEMP1, TEMP2 ) );

            // Choose random values for IL and IU, and ask for the
            // IL-th through IU-th eigenvalues.

            NTEST = 19;
            if ( N <= 1 ) {
               IL = 1;
               IU = N;
            } else {
               IL = 1 + ( N-1 )*INT( DLARND( 1, ISEED2 ) );
               IU = 1 + ( N-1 )*INT( DLARND( 1, ISEED2 ) );
               if ( IU < IL ) {
                  ITEMP = IU;
                  IU = IL;
                  IL = ITEMP;
               }
            }

            dstebz('I', 'E', N, VL, VU, IL, IU, ABSTOL, SD, SE, M2, NSPLIT, WA2, IWORK( 1 ), IWORK( N+1 ), WORK, IWORK( 2*N+1 ), IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'DSTEBZ(I)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               if ( IINFO < 0 ) {
                  return;
               } else {
                  RESULT[19] = ULPINV;
                  GO TO 280;
               }
            }

            // Determine the values VL and VU of the IL-th and IU-th
            // eigenvalues and ask for all eigenvalues in this range.

            if ( N > 0 ) {
               if ( IL != 1 ) {
                  VL = WA1( IL ) - max( HALF*( WA1( IL )-WA1( IL-1 ) ), ULP*ANORM, TWO*RTUNFL );
               } else {
                  VL = WA1( 1 ) - max( HALF*( WA1( N )-WA1( 1 ) ), ULP*ANORM, TWO*RTUNFL );
               }
               if ( IU != N ) {
                  VU = WA1( IU ) + max( HALF*( WA1( IU+1 )-WA1( IU ) ), ULP*ANORM, TWO*RTUNFL );
               } else {
                  VU = WA1( N ) + max( HALF*( WA1( N )-WA1( 1 ) ), ULP*ANORM, TWO*RTUNFL );
               }
            } else {
               VL = ZERO;
               VU = ONE;
            }

            dstebz('V', 'E', N, VL, VU, IL, IU, ABSTOL, SD, SE, M3, NSPLIT, WA3, IWORK( 1 ), IWORK( N+1 ), WORK, IWORK( 2*N+1 ), IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'DSTEBZ(V)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               if ( IINFO < 0 ) {
                  return;
               } else {
                  RESULT[19] = ULPINV;
                  GO TO 280;
               }
            }

            if ( M3 == 0 && N != 0 ) {
               RESULT[19] = ULPINV;
               GO TO 280;
            }

            // Do test 19

            TEMP1 = DSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL );
            TEMP2 = DSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL );
            if ( N > 0 ) {
               TEMP3 = max( ( WA1( N ) ).abs(), ( WA1( 1 ) ) ).abs();
            } else {
               TEMP3 = ZERO;
            }

            RESULT[19] = ( TEMP1+TEMP2 ) / max( UNFL, TEMP3*ULP );

            // Call DSTEIN to compute eigenvectors corresponding to
            // eigenvalues in WA1.  (First call DSTEBZ again, to make sure
            // it returns these eigenvalues in the correct order.)

            NTEST = 21;
            dstebz('A', 'B', N, VL, VU, IL, IU, ABSTOL, SD, SE, M, NSPLIT, WA1, IWORK( 1 ), IWORK( N+1 ), WORK, IWORK( 2*N+1 ), IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'DSTEBZ(A,B)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               if ( IINFO < 0 ) {
                  return;
               } else {
                  RESULT[20] = ULPINV;
                  RESULT[21] = ULPINV;
                  GO TO 280;
               }
            }

            dstein(N, SD, SE, M, WA1, IWORK( 1 ), IWORK( N+1 ), Z, LDU, WORK, IWORK( 2*N+1 ), IWORK( 3*N+1 ), IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'DSTEIN', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               if ( IINFO < 0 ) {
                  return;
               } else {
                  RESULT[20] = ULPINV;
                  RESULT[21] = ULPINV;
                  GO TO 280;
               }
            }

            // Do tests 20 and 21

            dstt21(N, 0, SD, SE, WA1, DUMMA, Z, LDU, WORK, RESULT( 20 ) );

            // Call DSTEDC(I) to compute D1 and Z, do tests.

            // Compute D1 and Z

            dcopy(N, SD, 1, D1, 1 );
            if (N > 0) dcopy( N-1, SE, 1, WORK, 1 );
            dlaset('Full', N, N, ZERO, ONE, Z, LDU );

            NTEST = 22;
            dstedc('I', N, D1, WORK, Z, LDU, WORK( N+1 ), LWEDC-N, IWORK, LIWEDC, IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'DSTEDC(I)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               if ( IINFO < 0 ) {
                  return;
               } else {
                  RESULT[22] = ULPINV;
                  GO TO 280;
               }
            }

            // Do Tests 22 and 23

            dstt21(N, 0, SD, SE, D1, DUMMA, Z, LDU, WORK, RESULT( 22 ) );

            // Call DSTEDC(V) to compute D1 and Z, do tests.

            // Compute D1 and Z

            dcopy(N, SD, 1, D1, 1 );
            if (N > 0) dcopy( N-1, SE, 1, WORK, 1 );
            dlaset('Full', N, N, ZERO, ONE, Z, LDU );

            NTEST = 24;
            dstedc('V', N, D1, WORK, Z, LDU, WORK( N+1 ), LWEDC-N, IWORK, LIWEDC, IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'DSTEDC(V)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               if ( IINFO < 0 ) {
                  return;
               } else {
                  RESULT[24] = ULPINV;
                  GO TO 280;
               }
            }

            // Do Tests 24 and 25

            dstt21(N, 0, SD, SE, D1, DUMMA, Z, LDU, WORK, RESULT( 24 ) );

            // Call DSTEDC(N) to compute D2, do tests.

            // Compute D2

            dcopy(N, SD, 1, D2, 1 );
            if (N > 0) dcopy( N-1, SE, 1, WORK, 1 );
            dlaset('Full', N, N, ZERO, ONE, Z, LDU );

            NTEST = 26;
            dstedc('N', N, D2, WORK, Z, LDU, WORK( N+1 ), LWEDC-N, IWORK, LIWEDC, IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'DSTEDC(N)', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               if ( IINFO < 0 ) {
                  return;
               } else {
                  RESULT[26] = ULPINV;
                  GO TO 280;
               }
            }

            // Do Test 26

            TEMP1 = ZERO;
            TEMP2 = ZERO;

            for (J = 1; J <= N; J++) { // 210
               TEMP1 = max( TEMP1, ( D1( J ) ).abs(), ( D2( J ) ) ).abs();
               TEMP2 = max( TEMP2, ABS( D1( J )-D2( J ) ) );
            } // 210

            RESULT[26] = TEMP2 / max( UNFL, ULP*max( TEMP1, TEMP2 ) );

            // Only test DSTEMR if IEEE compliant

            if ( ILAENV( 10, 'DSTEMR', 'VA', 1, 0, 0, 0 ) == 1 && ILAENV( 11, 'DSTEMR', 'VA', 1, 0, 0, 0 ) == 1 ) {

            // Call DSTEMR, do test 27 (relative eigenvalue accuracy)

               // If S is positive definite and diagonally dominant,
               // ask for all eigenvalues with high relative accuracy.

               VL = ZERO;
               VU = ZERO;
               IL = 0;
               IU = 0;
               if ( JTYPE == 21 && SREL ) {
                  NTEST = 27;
                  ABSTOL = UNFL + UNFL;
                  dstemr('V', 'A', N, SD, SE, VL, VU, IL, IU, M, WR, Z, LDU, N, IWORK( 1 ), TRYRAC, WORK, LWORK, IWORK( 2*N+1 ), LWORK-2*N, IINFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'DSTEMR(V,A,rel)', IINFO, N, JTYPE, IOLDSD;
                     INFO = ( IINFO ).abs();
                     if ( IINFO < 0 ) {
                        return;
                     } else {
                        RESULT[27] = ULPINV;
                        GO TO 270;
                     }
                  }

               // Do test 27

                  TEMP2 = TWO*( TWO*N-ONE )*ULP*( ONE+EIGHT*HALF**2 ) / ( ONE-HALF )**4;

                  TEMP1 = ZERO;
                  for (J = 1; J <= N; J++) { // 220
                     TEMP1 = max( TEMP1, ABS( D4( J )-WR( N-J+1 ) ) / ( ABSTOL+( D4( J ) ) ) ).abs();
                  } // 220

                  RESULT[27] = TEMP1 / TEMP2;

                  IL = 1 + ( N-1 )*INT( DLARND( 1, ISEED2 ) );
                  IU = 1 + ( N-1 )*INT( DLARND( 1, ISEED2 ) );
                  if ( IU < IL ) {
                     ITEMP = IU;
                     IU = IL;
                     IL = ITEMP;
                  }

                  if ( SRANGE ) {
                     NTEST = 28;
                     ABSTOL = UNFL + UNFL;
                     dstemr('V', 'I', N, SD, SE, VL, VU, IL, IU, M, WR, Z, LDU, N, IWORK( 1 ), TRYRAC, WORK, LWORK, IWORK( 2*N+1 ), LWORK-2*N, IINFO );

                     if ( IINFO != 0 ) {
                        WRITE( NOUNIT, FMT = 9999 )'DSTEMR(V,I,rel)', IINFO, N, JTYPE, IOLDSD;
                        INFO = ( IINFO ).abs();
                        if ( IINFO < 0 ) {
                           return;
                        } else {
                           RESULT[28] = ULPINV;
                           GO TO 270;
                        }
                     }

                  // Do test 28

                     TEMP2 = TWO*( TWO*N-ONE )*ULP* ( ONE+EIGHT*HALF**2 ) / ( ONE-HALF )**4;

                     TEMP1 = ZERO;
                     for (J = IL; J <= IU; J++) { // 230
                        TEMP1 = max( TEMP1, ABS( WR( J-IL+1 )-D4( N-J+ 1 ) ) / ( ABSTOL+( WR( J-IL+1 ) ) ) ).abs();
                     } // 230

                     RESULT[28] = TEMP1 / TEMP2;
                  } else {
                     RESULT[28] = ZERO;
                  }
               } else {
                  RESULT[27] = ZERO;
                  RESULT[28] = ZERO;
               }

            // Call DSTEMR(V,I) to compute D1 and Z, do tests.

            // Compute D1 and Z

               dcopy(N, SD, 1, D5, 1 );
               if (N > 0) dcopy( N-1, SE, 1, WORK, 1 );
               dlaset('Full', N, N, ZERO, ONE, Z, LDU );

               if ( SRANGE ) {
                  NTEST = 29;
                  IL = 1 + ( N-1 )*INT( DLARND( 1, ISEED2 ) );
                  IU = 1 + ( N-1 )*INT( DLARND( 1, ISEED2 ) );
                  if ( IU < IL ) {
                     ITEMP = IU;
                     IU = IL;
                     IL = ITEMP;
                  }
                  dstemr('V', 'I', N, D5, WORK, VL, VU, IL, IU, M, D1, Z, LDU, N, IWORK( 1 ), TRYRAC, WORK( N+1 ), LWORK-N, IWORK( 2*N+1 ), LIWORK-2*N, IINFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'DSTEMR(V,I)', IINFO, N, JTYPE, IOLDSD;
                     INFO = ( IINFO ).abs();
                     if ( IINFO < 0 ) {
                        return;
                     } else {
                        RESULT[29] = ULPINV;
                        GO TO 280;
                     }
                  }

            // Do Tests 29 and 30

                  dstt22(N, M, 0, SD, SE, D1, DUMMA, Z, LDU, WORK, M, RESULT( 29 ) );

            // Call DSTEMR to compute D2, do tests.

            // Compute D2

                  dcopy(N, SD, 1, D5, 1 );
                  if (N > 0) dcopy( N-1, SE, 1, WORK, 1 );

                  NTEST = 31;
                  dstemr('N', 'I', N, D5, WORK, VL, VU, IL, IU, M, D2, Z, LDU, N, IWORK( 1 ), TRYRAC, WORK( N+1 ), LWORK-N, IWORK( 2*N+1 ), LIWORK-2*N, IINFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'DSTEMR(N,I)', IINFO, N, JTYPE, IOLDSD;
                     INFO = ( IINFO ).abs();
                     if ( IINFO < 0 ) {
                        return;
                     } else {
                        RESULT[31] = ULPINV;
                        GO TO 280;
                     }
                  }

            // Do Test 31

                  TEMP1 = ZERO;
                  TEMP2 = ZERO;

                  for (J = 1; J <= IU - IL + 1; J++) { // 240
                     TEMP1 = max( TEMP1, ( D1( J ) ).abs(), ( D2( J ) ) ).abs();
                     TEMP2 = max( TEMP2, ABS( D1( J )-D2( J ) ) );
                  } // 240

                  RESULT[31] = TEMP2 / max( UNFL, ULP*max( TEMP1, TEMP2 ) );

            // Call DSTEMR(V,V) to compute D1 and Z, do tests.

            // Compute D1 and Z

                  dcopy(N, SD, 1, D5, 1 );
                  if (N > 0) dcopy( N-1, SE, 1, WORK, 1 );
                  dlaset('Full', N, N, ZERO, ONE, Z, LDU );

                  NTEST = 32;

                  if ( N > 0 ) {
                     if ( IL != 1 ) {
                        VL = D2( IL ) - max( HALF* ( D2( IL )-D2( IL-1 ) ), ULP*ANORM, TWO*RTUNFL );
                     } else {
                        VL = D2( 1 ) - max( HALF*( D2( N )-D2( 1 ) ), ULP*ANORM, TWO*RTUNFL );
                     }
                     if ( IU != N ) {
                        VU = D2( IU ) + max( HALF* ( D2( IU+1 )-D2( IU ) ), ULP*ANORM, TWO*RTUNFL );
                     } else {
                        VU = D2( N ) + max( HALF*( D2( N )-D2( 1 ) ), ULP*ANORM, TWO*RTUNFL );
                     }
                  } else {
                     VL = ZERO;
                     VU = ONE;
                  }

                  dstemr('V', 'V', N, D5, WORK, VL, VU, IL, IU, M, D1, Z, LDU, N, IWORK( 1 ), TRYRAC, WORK( N+1 ), LWORK-N, IWORK( 2*N+1 ), LIWORK-2*N, IINFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'DSTEMR(V,V)', IINFO, N, JTYPE, IOLDSD;
                     INFO = ( IINFO ).abs();
                     if ( IINFO < 0 ) {
                        return;
                     } else {
                        RESULT[32] = ULPINV;
                        GO TO 280;
                     }
                  }

            // Do Tests 32 and 33

                  dstt22(N, M, 0, SD, SE, D1, DUMMA, Z, LDU, WORK, M, RESULT( 32 ) );

            // Call DSTEMR to compute D2, do tests.

            // Compute D2

                  dcopy(N, SD, 1, D5, 1 );
                  if (N > 0) dcopy( N-1, SE, 1, WORK, 1 );

                  NTEST = 34;
                  dstemr('N', 'V', N, D5, WORK, VL, VU, IL, IU, M, D2, Z, LDU, N, IWORK( 1 ), TRYRAC, WORK( N+1 ), LWORK-N, IWORK( 2*N+1 ), LIWORK-2*N, IINFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'DSTEMR(N,V)', IINFO, N, JTYPE, IOLDSD;
                     INFO = ( IINFO ).abs();
                     if ( IINFO < 0 ) {
                        return;
                     } else {
                        RESULT[34] = ULPINV;
                        GO TO 280;
                     }
                  }

            // Do Test 34

                  TEMP1 = ZERO;
                  TEMP2 = ZERO;

                  for (J = 1; J <= IU - IL + 1; J++) { // 250
                     TEMP1 = max( TEMP1, ( D1( J ) ).abs(), ( D2( J ) ) ).abs();
                     TEMP2 = max( TEMP2, ABS( D1( J )-D2( J ) ) );
                  } // 250

                  RESULT[34] = TEMP2 / max( UNFL, ULP*max( TEMP1, TEMP2 ) );
               } else {
                  RESULT[29] = ZERO;
                  RESULT[30] = ZERO;
                  RESULT[31] = ZERO;
                  RESULT[32] = ZERO;
                  RESULT[33] = ZERO;
                  RESULT[34] = ZERO;
               }

            // Call DSTEMR(V,A) to compute D1 and Z, do tests.

            // Compute D1 and Z

               dcopy(N, SD, 1, D5, 1 );
               if (N > 0) dcopy( N-1, SE, 1, WORK, 1 );

               NTEST = 35;

               dstemr('V', 'A', N, D5, WORK, VL, VU, IL, IU, M, D1, Z, LDU, N, IWORK( 1 ), TRYRAC, WORK( N+1 ), LWORK-N, IWORK( 2*N+1 ), LIWORK-2*N, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'DSTEMR(V,A)', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[35] = ULPINV;
                     GO TO 280;
                  }
               }

            // Do Tests 35 and 36

               dstt22(N, M, 0, SD, SE, D1, DUMMA, Z, LDU, WORK, M, RESULT( 35 ) );

            // Call DSTEMR to compute D2, do tests.

            // Compute D2

               dcopy(N, SD, 1, D5, 1 );
               if (N > 0) dcopy( N-1, SE, 1, WORK, 1 );

               NTEST = 37;
               dstemr('N', 'A', N, D5, WORK, VL, VU, IL, IU, M, D2, Z, LDU, N, IWORK( 1 ), TRYRAC, WORK( N+1 ), LWORK-N, IWORK( 2*N+1 ), LIWORK-2*N, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'DSTEMR(N,A)', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[37] = ULPINV;
                     GO TO 280;
                  }
               }

            // Do Test 37

               TEMP1 = ZERO;
               TEMP2 = ZERO;

               for (J = 1; J <= N; J++) { // 260
                  TEMP1 = max( TEMP1, ( D1( J ) ).abs(), ( D2( J ) ) ).abs();
                  TEMP2 = max( TEMP2, ABS( D1( J )-D2( J ) ) );
               } // 260

               RESULT[37] = TEMP2 / max( UNFL, ULP*max( TEMP1, TEMP2 ) );
            }
            } // 270
            } // 280
            NTESTT = NTESTT + NTEST;

            // End of Loop -- Check for RESULT(j) > THRESH

            // Print out tests which fail.

            for (JR = 1; JR <= NTEST; JR++) { // 290
               if ( RESULT( JR ) >= THRESH ) {

                  // If this is the first test to fail,
                  // print a header to the data file.

                  if ( NERRS == 0 ) {
                     WRITE( NOUNIT, FMT = 9998 )'DST';
                     WRITE( NOUNIT, FMT = 9997 );
                     WRITE( NOUNIT, FMT = 9996 );
                     WRITE( NOUNIT, FMT = 9995 )'Symmetric';
                     WRITE( NOUNIT, FMT = 9994 );

                     // Tests performed

                     WRITE( NOUNIT, FMT = 9988 );
                  }
                  NERRS = NERRS + 1;
                  WRITE( NOUNIT, FMT = 9990 )N, IOLDSD, JTYPE, JR, RESULT( JR );
               }
            } // 290
         } // 300
      } // 310

      // Summary

      dlasum('DST', NOUNIT, NERRS, NTESTT );
      return;

 9999 FORMAT( ' DCHKST2STG: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' );

 9998 FORMAT( / 1X, A3, ' -- Real Symmetric eigenvalue problem' );
 9997 FORMAT( ' Matrix types (see DCHKST2STG for details): ' );

 9996 FORMAT( / ' Special Matrices:', / '  1=Zero matrix.                        ', '  5=Diagonal: clustered entries.', / '  2=Identity matrix.                    ', '  6=Diagonal: large, evenly spaced.', / '  3=Diagonal: evenly spaced entries.    ', '  7=Diagonal: small, evenly spaced.', / '  4=Diagonal: geometr. spaced entries.' );
 9995 FORMAT( ' Dense ', A, ' Matrices:', / '  8=Evenly spaced eigenvals.            ', ' 12=Small, evenly spaced eigenvals.', / '  9=Geometrically spaced eigenvals.     ', ' 13=Matrix with random O(1) entries.', / ' 10=Clustered eigenvalues.              ', ' 14=Matrix with large random entries.', / ' 11=Large, evenly spaced eigenvals.     ', ' 15=Matrix with small random entries.' );
 9994 FORMAT( ' 16=Positive definite, evenly spaced eigenvalues', / ' 17=Positive definite, geometrically spaced eigenvlaues', / ' 18=Positive definite, clustered eigenvalues', / ' 19=Positive definite, small evenly spaced eigenvalues', / ' 20=Positive definite, large evenly spaced eigenvalues', / ' 21=Diagonally dominant tridiagonal, geometrically', ' spaced eigenvalues' );

 9990 FORMAT( ' N=', I5, ', seed=', 4( I4, ',' ), ' type ', I2, ', test(', I2, ')=', G10.3 );

 9988 FORMAT( / 'Test performed:  see DCHKST2STG for details.', / );
      }