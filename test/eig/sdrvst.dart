      void sdrvst(NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH, NOUNIT, A, LDA, D1, D2, D3, D4, EVEIGS, WA1, WA2, WA3, U, LDU, V, TAU, Z, WORK, LWORK, IWORK, LIWORK, RESULT, INFO ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDU, LIWORK, LWORK, NOUNIT, NSIZES, NTYPES;
      double               THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                ISEED( 4 ), IWORK( * ), NN( * );
      double               A( LDA, * ), D1( * ), D2( * ), D3( * ), D4( * ), EVEIGS( * ), RESULT( * ), TAU( * ), U( LDU, * ), V( LDU, * ), WA1( * ), WA2( * ), WA3( * ), WORK( * ), Z( LDU, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double               ZERO, ONE, TWO, TEN;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0, TEN = 10.0 ;
      double               HALF;
      const              HALF = 0.5 ;
      int                MAXTYP;
      const              MAXTYP = 18 ;
      // ..
      // .. Local Scalars ..
      bool               BADNN;
      String             UPLO;
      int                I, IDIAG, IHBW, IINFO, IL, IMODE, INDX, IROW, ITEMP, ITYPE, IU, IUPLO, J, J1, J2, JCOL, JSIZE, JTYPE, KD, LGN, LIWEDC, LWEDC, M, M2, M3, MTYPES, N, NERRS, NMATS, NMAX, NTEST, NTESTT;
      double               ABSTOL, ANINV, ANORM, COND, OVFL, RTOVFL, RTUNFL, TEMP1, TEMP2, TEMP3, ULP, ULPINV, UNFL, VL, VU;
      // ..
      // .. Local Arrays ..
      int                IDUMMA( 1 ), IOLDSD( 4 ), ISEED2( 4 ), ISEED3( 4 ), KMAGN( MAXTYP ), KMODE( MAXTYP ), KTYPE( MAXTYP );
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH, SLARND, SSXT1;
      // EXTERNAL SLAMCH, SLARND, SSXT1
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALASVM, SLACPY, SLAFTS, SLASET, SLATMR, SLATMS, SSBEV, SSBEVD, SSBEVX, SSPEV, SSPEVD, SSPEVX, SSTEV, SSTEVD, SSTEVR, SSTEVX, SSTT21, SSTT22, SSYEV, SSYEVD, SSYEVR, SSYEVX, SSYT21, SSYT22, XERBLA
      // ..
      // .. Scalars in Common ..
      String             SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, INT, LOG, MAX, MIN, REAL, SQRT
      // ..
      // .. Data statements ..
      const KTYPE = [ 1, 2, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 8, 8, 8, 9, 9, 9,];
      const KMAGN = [ 1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3, 1, 2, 3 ];
      const KMODE = [ 0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0, 4, 4, 4 ];
      // ..
      // .. Executable Statements ..

      // Keep ftrnchek happy

      VL = ZERO;
      VU = ZERO;

      // 1)      Check for errors

      NTESTT = 0;
      INFO = 0;

      BADNN = false;
      NMAX = 1;
      for (J = 1; J <= NSIZES; J++) { // 10
         NMAX = max( NMAX, NN( J ) );
         if( NN( J ) < 0 ) BADNN = true;
      } // 10

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
         INFO = -16;
      } else if ( 2*max( 2, NMAX )**2 > LWORK ) {
         INFO = -21;
      }

      if ( INFO != 0 ) {
         xerbla('SDRVST', -INFO );
         return;
      }

      // Quick return if nothing to do

      if (NSIZES == 0 || NTYPES == 0) return;

      // More Important constants

      UNFL = SLAMCH( 'Safe minimum' );
      OVFL = SLAMCH( 'Overflow' );
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' );
      ULPINV = ONE / ULP;
      RTUNFL = sqrt( UNFL );
      RTOVFL = sqrt( OVFL );

      // Loop over sizes, types

      for (I = 1; I <= 4; I++) { // 20
         ISEED2[I] = ISEED( I );
         ISEED3[I] = ISEED( I );
      } // 20

      NERRS = 0;
      NMATS = 0;


      for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) { // 1740
         N = NN( JSIZE );
         if ( N > 0 ) {
            LGN = INT( LOG( double( N ) ) / LOG( TWO ) );
            if (2**LGN < N) LGN = LGN + 1;
            IF( 2**LGN < N ) LGN = LGN + 1;
            LWEDC = 1 + 4*N + 2*N*LGN + 4*N**2;
            // LIWEDC = 6 + 6*N + 5*N*LGN
            LIWEDC = 3 + 5*N;
         } else {
            LWEDC = 9;
            // LIWEDC = 12
            LIWEDC = 8;
         }
         ANINV = ONE / REAL( max( 1, N ) );

         if ( NSIZES != 1 ) {
            MTYPES = min( MAXTYP, NTYPES );
         } else {
            MTYPES = min( MAXTYP+1, NTYPES );
         }

         for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) { // 1730

            if( !DOTYPE( JTYPE ) ) GO TO 1730;
            NMATS = NMATS + 1;
            NTEST = 0;

            for (J = 1; J <= 4; J++) { // 30
               IOLDSD[J] = ISEED( J );
            } // 30

            // 2)      Compute "A"

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
            // =9                      band symmetric, w/ eigenvalues

            if (MTYPES > MAXTYP) GO TO 110;

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

            slaset('Full', LDA, N, ZERO, ZERO, A, LDA );
            IINFO = 0;
            COND = ULPINV;

            // Special Matrices -- Identity & Jordan block

                    // Zero

            if ( ITYPE == 1 ) {
               IINFO = 0;

            } else if ( ITYPE == 2 ) {

               // Identity

               for (JCOL = 1; JCOL <= N; JCOL++) { // 80
                  A[JCOL, JCOL] = ANORM;
               } // 80

            } else if ( ITYPE == 4 ) {

               // Diagonal Matrix, [Eigen]values Specified

               slatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, 0, 0, 'N', A, LDA, WORK( N+1 ), IINFO );

            } else if ( ITYPE == 5 ) {

               // Symmetric, eigenvalues specified

               slatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, N, N, 'N', A, LDA, WORK( N+1 ), IINFO );

            } else if ( ITYPE == 7 ) {

               // Diagonal, random eigenvalues

               IDUMMA[1] = 1;
               slatmr(N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, 0, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else if ( ITYPE == 8 ) {

               // Symmetric, random eigenvalues

               IDUMMA[1] = 1;
               slatmr(N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else if ( ITYPE == 9 ) {

               // Symmetric banded, eigenvalues specified

               IHBW = INT( ( N-1 )*SLARND( 1, ISEED3 ) );
               slatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, IHBW, IHBW, 'Z', U, LDU, WORK( N+1 ), IINFO );

               // Store as dense matrix for most routines.

               slaset('Full', LDA, N, ZERO, ZERO, A, LDA );
               for (IDIAG = -IHBW; IDIAG <= IHBW; IDIAG++) { // 100
                  IROW = IHBW - IDIAG + 1;
                  J1 = max( 1, IDIAG+1 );
                  J2 = min( N, N+IDIAG );
                  for (J = J1; J <= J2; J++) { // 90
                     I = J - IDIAG;
                     A[I, J] = U( IROW, J );
                  } // 90
               } // 100
            } else {
               IINFO = 1;
            }

            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'Generator', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               return;
            }

            } // 110

            ABSTOL = UNFL + UNFL;
            if ( N <= 1 ) {
               IL = 1;
               IU = N;
            } else {
               IL = 1 + INT( ( N-1 )*SLARND( 1, ISEED2 ) );
               IU = 1 + INT( ( N-1 )*SLARND( 1, ISEED2 ) );
               if ( IL > IU ) {
                  ITEMP = IL;
                  IL = IU;
                  IU = ITEMP;
               }
            }

            // 3)      If matrix is tridiagonal, call SSTEV and SSTEVX.

            if ( JTYPE <= 7 ) {
               NTEST = 1;
               for (I = 1; I <= N; I++) { // 120
                  D1[I] = double( A( I, I ) );
               } // 120
               for (I = 1; I <= N - 1; I++) { // 130
                  D2[I] = double( A( I+1, I ) );
               } // 130
               SRNAMT = 'SSTEV';
               sstev('V', N, D1, D2, Z, LDU, WORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEV(V)', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[1] = ULPINV;
                     RESULT[2] = ULPINV;
                     RESULT[3] = ULPINV;
                     GO TO 180;
                  }
               }

               // Do tests 1 and 2.

               for (I = 1; I <= N; I++) { // 140
                  D3[I] = double( A( I, I ) );
               } // 140
               for (I = 1; I <= N - 1; I++) { // 150
                  D4[I] = double( A( I+1, I ) );
               } // 150
               sstt21(N, 0, D3, D4, D1, D2, Z, LDU, WORK, RESULT( 1 ) );

               NTEST = 3;
               for (I = 1; I <= N - 1; I++) { // 160
                  D4[I] = double( A( I+1, I ) );
               } // 160
               SRNAMT = 'SSTEV';
               sstev('N', N, D3, D4, Z, LDU, WORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEV(N)', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[3] = ULPINV;
                     GO TO 180;
                  }
               }

               // Do test 3.

               TEMP1 = ZERO;
               TEMP2 = ZERO;
               for (J = 1; J <= N; J++) { // 170
                  TEMP1 = max( TEMP1, ( D1( J ) ).abs(), ( D3( J ) ) ).abs();
                  TEMP2 = max( TEMP2, ABS( D1( J )-D3( J ) ) );
               } // 170
               RESULT[3] = TEMP2 / max( UNFL, ULP*max( TEMP1, TEMP2 ) );

               } // 180

               NTEST = 4;
               for (I = 1; I <= N; I++) { // 190
                  EVEIGS[I] = D3( I );
                  D1[I] = double( A( I, I ) );
               } // 190
               for (I = 1; I <= N - 1; I++) { // 200
                  D2[I] = double( A( I+1, I ) );
               } // 200
               SRNAMT = 'SSTEVX';
               sstevx('V', 'A', N, D1, D2, VL, VU, IL, IU, ABSTOL, M, WA1, Z, LDU, WORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEVX(V,A)', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[4] = ULPINV;
                     RESULT[5] = ULPINV;
                     RESULT[6] = ULPINV;
                     GO TO 250;
                  }
               }
               if ( N > 0 ) {
                  TEMP3 = max( ( WA1( 1 ) ).abs(), ( WA1( N ) ) ).abs();
               } else {
                  TEMP3 = ZERO;
               }

               // Do tests 4 and 5.

               for (I = 1; I <= N; I++) { // 210
                  D3[I] = double( A( I, I ) );
               } // 210
               for (I = 1; I <= N - 1; I++) { // 220
                  D4[I] = double( A( I+1, I ) );
               } // 220
               sstt21(N, 0, D3, D4, WA1, D2, Z, LDU, WORK, RESULT( 4 ) );

               NTEST = 6;
               for (I = 1; I <= N - 1; I++) { // 230
                  D4[I] = double( A( I+1, I ) );
               } // 230
               SRNAMT = 'SSTEVX';
               sstevx('N', 'A', N, D3, D4, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEVX(N,A)', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[6] = ULPINV;
                     GO TO 250;
                  }
               }

               // Do test 6.

               TEMP1 = ZERO;
               TEMP2 = ZERO;
               for (J = 1; J <= N; J++) { // 240
                  TEMP1 = max( TEMP1, ( WA2( J ) ).abs(), ( EVEIGS( J ) ) ).abs();
                  TEMP2 = max( TEMP2, ABS( WA2( J )-EVEIGS( J ) ) );
               } // 240
               RESULT[6] = TEMP2 / max( UNFL, ULP*max( TEMP1, TEMP2 ) );

               } // 250

               NTEST = 7;
               for (I = 1; I <= N; I++) { // 260
                  D1[I] = double( A( I, I ) );
               } // 260
               for (I = 1; I <= N - 1; I++) { // 270
                  D2[I] = double( A( I+1, I ) );
               } // 270
               SRNAMT = 'SSTEVR';
               sstevr('V', 'A', N, D1, D2, VL, VU, IL, IU, ABSTOL, M, WA1, Z, LDU, IWORK, WORK, LWORK, IWORK(2*N+1), LIWORK-2*N, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEVR(V,A)', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[7] = ULPINV;
                     RESULT[8] = ULPINV;
                     GO TO 320;
                  }
               }
               if ( N > 0 ) {
                  TEMP3 = max( ( WA1( 1 ) ).abs(), ( WA1( N ) ) ).abs();
               } else {
                  TEMP3 = ZERO;
               }

               // Do tests 7 and 8.

               for (I = 1; I <= N; I++) { // 280
                  D3[I] = double( A( I, I ) );
               } // 280
               for (I = 1; I <= N - 1; I++) { // 290
                  D4[I] = double( A( I+1, I ) );
               } // 290
               sstt21(N, 0, D3, D4, WA1, D2, Z, LDU, WORK, RESULT( 7 ) );

               NTEST = 9;
               for (I = 1; I <= N - 1; I++) { // 300
                  D4[I] = double( A( I+1, I ) );
               } // 300
               SRNAMT = 'SSTEVR';
               sstevr('N', 'A', N, D3, D4, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, IWORK, WORK, LWORK, IWORK(2*N+1), LIWORK-2*N, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEVR(N,A)', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[9] = ULPINV;
                     GO TO 320;
                  }
               }

               // Do test 9.

               TEMP1 = ZERO;
               TEMP2 = ZERO;
               for (J = 1; J <= N; J++) { // 310
                  TEMP1 = max( TEMP1, ( WA2( J ) ).abs(), ( EVEIGS( J ) ) ).abs();
                  TEMP2 = max( TEMP2, ABS( WA2( J )-EVEIGS( J ) ) );
               } // 310
               RESULT[9] = TEMP2 / max( UNFL, ULP*max( TEMP1, TEMP2 ) );

               } // 320


               NTEST = 10;
               for (I = 1; I <= N; I++) { // 330
                  D1[I] = double( A( I, I ) );
               } // 330
               for (I = 1; I <= N - 1; I++) { // 340
                  D2[I] = double( A( I+1, I ) );
               } // 340
               SRNAMT = 'SSTEVX';
               sstevx('V', 'I', N, D1, D2, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEVX(V,I)', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[10] = ULPINV;
                     RESULT[11] = ULPINV;
                     RESULT[12] = ULPINV;
                     GO TO 380;
                  }
               }

               // Do tests 10 and 11.

               for (I = 1; I <= N; I++) { // 350
                  D3[I] = double( A( I, I ) );
               } // 350
               for (I = 1; I <= N - 1; I++) { // 360
                  D4[I] = double( A( I+1, I ) );
               } // 360
               sstt22(N, M2, 0, D3, D4, WA2, D2, Z, LDU, WORK, max( 1, M2 ), RESULT( 10 ) );


               NTEST = 12;
               for (I = 1; I <= N - 1; I++) { // 370
                  D4[I] = double( A( I+1, I ) );
               } // 370
               SRNAMT = 'SSTEVX';
               sstevx('N', 'I', N, D3, D4, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, WORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEVX(N,I)', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[12] = ULPINV;
                     GO TO 380;
                  }
               }

               // Do test 12.

               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL );
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL );
               RESULT[12] = ( TEMP1+TEMP2 ) / max( UNFL, ULP*TEMP3 );

               } // 380

               NTEST = 12;
               if ( N > 0 ) {
                  if ( IL != 1 ) {
                     VL = WA1( IL ) - max( HALF* ( WA1( IL )-WA1( IL-1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL );
                  } else {
                     VL = WA1( 1 ) - max( HALF*( WA1( N )-WA1( 1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL );
                  }
                  if ( IU != N ) {
                     VU = WA1( IU ) + max( HALF* ( WA1( IU+1 )-WA1( IU ) ), TEN*ULP*TEMP3, TEN*RTUNFL );
                  } else {
                     VU = WA1( N ) + max( HALF*( WA1( N )-WA1( 1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL );
                  }
               } else {
                  VL = ZERO;
                  VU = ONE;
               }

               for (I = 1; I <= N; I++) { // 390
                  D1[I] = double( A( I, I ) );
               } // 390
               for (I = 1; I <= N - 1; I++) { // 400
                  D2[I] = double( A( I+1, I ) );
               } // 400
               SRNAMT = 'SSTEVX';
               sstevx('V', 'V', N, D1, D2, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEVX(V,V)', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[13] = ULPINV;
                     RESULT[14] = ULPINV;
                     RESULT[15] = ULPINV;
                     GO TO 440;
                  }
               }

               if ( M2 == 0 && N > 0 ) {
                  RESULT[13] = ULPINV;
                  RESULT[14] = ULPINV;
                  RESULT[15] = ULPINV;
                  GO TO 440;
               }

               // Do tests 13 and 14.

               for (I = 1; I <= N; I++) { // 410
                  D3[I] = double( A( I, I ) );
               } // 410
               for (I = 1; I <= N - 1; I++) { // 420
                  D4[I] = double( A( I+1, I ) );
               } // 420
               sstt22(N, M2, 0, D3, D4, WA2, D2, Z, LDU, WORK, max( 1, M2 ), RESULT( 13 ) );

               NTEST = 15;
               for (I = 1; I <= N - 1; I++) { // 430
                  D4[I] = double( A( I+1, I ) );
               } // 430
               SRNAMT = 'SSTEVX';
               sstevx('N', 'V', N, D3, D4, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, WORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEVX(N,V)', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[15] = ULPINV;
                     GO TO 440;
                  }
               }

               // Do test 15.

               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL );
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL );
               RESULT[15] = ( TEMP1+TEMP2 ) / max( UNFL, TEMP3*ULP );

               } // 440

               NTEST = 16;
               for (I = 1; I <= N; I++) { // 450
                  D1[I] = double( A( I, I ) );
               } // 450
               for (I = 1; I <= N - 1; I++) { // 460
                  D2[I] = double( A( I+1, I ) );
               } // 460
               SRNAMT = 'SSTEVD';
               sstevd('V', N, D1, D2, Z, LDU, WORK, LWEDC, IWORK, LIWEDC, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEVD(V)', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[16] = ULPINV;
                     RESULT[17] = ULPINV;
                     RESULT[18] = ULPINV;
                     GO TO 510;
                  }
               }

               // Do tests 16 and 17.

               for (I = 1; I <= N; I++) { // 470
                  D3[I] = double( A( I, I ) );
               } // 470
               for (I = 1; I <= N - 1; I++) { // 480
                  D4[I] = double( A( I+1, I ) );
               } // 480
               sstt21(N, 0, D3, D4, D1, D2, Z, LDU, WORK, RESULT( 16 ) );

               NTEST = 18;
               for (I = 1; I <= N - 1; I++) { // 490
                  D4[I] = double( A( I+1, I ) );
               } // 490
               SRNAMT = 'SSTEVD';
               sstevd('N', N, D3, D4, Z, LDU, WORK, LWEDC, IWORK, LIWEDC, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEVD(N)', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[18] = ULPINV;
                     GO TO 510;
                  }
               }

               // Do test 18.

               TEMP1 = ZERO;
               TEMP2 = ZERO;
               for (J = 1; J <= N; J++) { // 500
                  TEMP1 = max( TEMP1, ( EVEIGS( J ) ).abs(), ( D3( J ) ) ).abs();
                  TEMP2 = max( TEMP2, ABS( EVEIGS( J )-D3( J ) ) );
               } // 500
               RESULT[18] = TEMP2 / max( UNFL, ULP*max( TEMP1, TEMP2 ) );

               } // 510

               NTEST = 19;
               for (I = 1; I <= N; I++) { // 520
                  D1[I] = double( A( I, I ) );
               } // 520
               for (I = 1; I <= N - 1; I++) { // 530
                  D2[I] = double( A( I+1, I ) );
               } // 530
               SRNAMT = 'SSTEVR';
               sstevr('V', 'I', N, D1, D2, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, IWORK, WORK, LWORK, IWORK(2*N+1), LIWORK-2*N, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEVR(V,I)', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[19] = ULPINV;
                     RESULT[20] = ULPINV;
                     RESULT[21] = ULPINV;
                     GO TO 570;
                  }
               }

               // DO tests 19 and 20.

               for (I = 1; I <= N; I++) { // 540
                  D3[I] = double( A( I, I ) );
               } // 540
               for (I = 1; I <= N - 1; I++) { // 550
                  D4[I] = double( A( I+1, I ) );
               } // 550
               sstt22(N, M2, 0, D3, D4, WA2, D2, Z, LDU, WORK, max( 1, M2 ), RESULT( 19 ) );


               NTEST = 21;
               for (I = 1; I <= N - 1; I++) { // 560
                  D4[I] = double( A( I+1, I ) );
               } // 560
               SRNAMT = 'SSTEVR';
               sstevr('N', 'I', N, D3, D4, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, IWORK, WORK, LWORK, IWORK(2*N+1), LIWORK-2*N, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEVR(N,I)', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[21] = ULPINV;
                     GO TO 570;
                  }
               }

               // Do test 21.

               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL );
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL );
               RESULT[21] = ( TEMP1+TEMP2 ) / max( UNFL, ULP*TEMP3 );

               } // 570

               NTEST = 21;
               if ( N > 0 ) {
                  if ( IL != 1 ) {
                     VL = WA1( IL ) - max( HALF* ( WA1( IL )-WA1( IL-1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL );
                  } else {
                     VL = WA1( 1 ) - max( HALF*( WA1( N )-WA1( 1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL );
                  }
                  if ( IU != N ) {
                     VU = WA1( IU ) + max( HALF* ( WA1( IU+1 )-WA1( IU ) ), TEN*ULP*TEMP3, TEN*RTUNFL );
                  } else {
                     VU = WA1( N ) + max( HALF*( WA1( N )-WA1( 1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL );
                  }
               } else {
                  VL = ZERO;
                  VU = ONE;
               }

               for (I = 1; I <= N; I++) { // 580
                  D1[I] = double( A( I, I ) );
               } // 580
               for (I = 1; I <= N - 1; I++) { // 590
                  D2[I] = double( A( I+1, I ) );
               } // 590
               SRNAMT = 'SSTEVR';
               sstevr('V', 'V', N, D1, D2, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, IWORK, WORK, LWORK, IWORK(2*N+1), LIWORK-2*N, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEVR(V,V)', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[22] = ULPINV;
                     RESULT[23] = ULPINV;
                     RESULT[24] = ULPINV;
                     GO TO 630;
                  }
               }

               if ( M2 == 0 && N > 0 ) {
                  RESULT[22] = ULPINV;
                  RESULT[23] = ULPINV;
                  RESULT[24] = ULPINV;
                  GO TO 630;
               }

               // Do tests 22 and 23.

               for (I = 1; I <= N; I++) { // 600
                  D3[I] = double( A( I, I ) );
               } // 600
               for (I = 1; I <= N - 1; I++) { // 610
                  D4[I] = double( A( I+1, I ) );
               } // 610
               sstt22(N, M2, 0, D3, D4, WA2, D2, Z, LDU, WORK, max( 1, M2 ), RESULT( 22 ) );

               NTEST = 24;
               for (I = 1; I <= N - 1; I++) { // 620
                  D4[I] = double( A( I+1, I ) );
               } // 620
               SRNAMT = 'SSTEVR';
               sstevr('N', 'V', N, D3, D4, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, IWORK, WORK, LWORK, IWORK(2*N+1), LIWORK-2*N, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSTEVR(N,V)', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[24] = ULPINV;
                     GO TO 630;
                  }
               }

               // Do test 24.

               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL );
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL );
               RESULT[24] = ( TEMP1+TEMP2 ) / max( UNFL, TEMP3*ULP );

               } // 630



            } else {

               for (I = 1; I <= 24; I++) { // 640
                  RESULT[I] = ZERO;
               } // 640
               NTEST = 24;
            }

            // Perform remaining tests storing upper or lower triangular
            // part of matrix.

            for (IUPLO = 0; IUPLO <= 1; IUPLO++) { // 1720
               if ( IUPLO == 0 ) {
                  UPLO = 'L';
               } else {
                  UPLO = 'U';
               }

               // 4)      Call SSYEV and SSYEVX.

               slacpy(' ', N, N, A, LDA, V, LDU );

               NTEST = NTEST + 1;
               SRNAMT = 'SSYEV';
               ssyev('V', UPLO, N, A, LDU, D1, WORK, LWORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSYEV(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     RESULT[NTEST+1] = ULPINV;
                     RESULT[NTEST+2] = ULPINV;
                     GO TO 660;
                  }
               }

               // Do tests 25 and 26 (or +54)

               ssyt21(1, UPLO, N, 0, V, LDU, D1, D2, A, LDU, Z, LDU, TAU, WORK, RESULT( NTEST ) );

               slacpy(' ', N, N, V, LDU, A, LDA );

               NTEST = NTEST + 2;
               SRNAMT = 'SSYEV';
               ssyev('N', UPLO, N, A, LDU, D3, WORK, LWORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSYEV(N,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     GO TO 660;
                  }
               }

               // Do test 27 (or +54)

               TEMP1 = ZERO;
               TEMP2 = ZERO;
               for (J = 1; J <= N; J++) { // 650
                  TEMP1 = max( TEMP1, ( D1( J ) ).abs(), ( D3( J ) ) ).abs();
                  TEMP2 = max( TEMP2, ABS( D1( J )-D3( J ) ) );
               } // 650
               RESULT[NTEST] = TEMP2 / max( UNFL, ULP*max( TEMP1, TEMP2 ) );

               } // 660
               slacpy(' ', N, N, V, LDU, A, LDA );

               NTEST = NTEST + 1;

               if ( N > 0 ) {
                  TEMP3 = max( ( D1( 1 ) ).abs(), ( D1( N ) ) ).abs();
                  if ( IL != 1 ) {
                     VL = D1( IL ) - max( HALF*( D1( IL )-D1( IL-1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL );
                  } else if ( N > 0 ) {
                     VL = D1( 1 ) - max( HALF*( D1( N )-D1( 1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL );
                  }
                  if ( IU != N ) {
                     VU = D1( IU ) + max( HALF*( D1( IU+1 )-D1( IU ) ), TEN*ULP*TEMP3, TEN*RTUNFL );
                  } else if ( N > 0 ) {
                     VU = D1( N ) + max( HALF*( D1( N )-D1( 1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL );
                  }
               } else {
                  TEMP3 = ZERO;
                  VL = ZERO;
                  VU = ONE;
               }

               SRNAMT = 'SSYEVX';
               ssyevx('V', 'A', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M, WA1, Z, LDU, WORK, LWORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSYEVX(V,A,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     RESULT[NTEST+1] = ULPINV;
                     RESULT[NTEST+2] = ULPINV;
                     GO TO 680;
                  }
               }

               // Do tests 28 and 29 (or +54)

               slacpy(' ', N, N, V, LDU, A, LDA );

               ssyt21(1, UPLO, N, 0, A, LDU, D1, D2, Z, LDU, V, LDU, TAU, WORK, RESULT( NTEST ) );

               NTEST = NTEST + 2;
               SRNAMT = 'SSYEVX';
               ssyevx('N', 'A', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK, LWORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSYEVX(N,A,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     GO TO 680;
                  }
               }

               // Do test 30 (or +54)

               TEMP1 = ZERO;
               TEMP2 = ZERO;
               for (J = 1; J <= N; J++) { // 670
                  TEMP1 = max( TEMP1, ( WA1( J ) ).abs(), ( WA2( J ) ) ).abs();
                  TEMP2 = max( TEMP2, ABS( WA1( J )-WA2( J ) ) );
               } // 670
               RESULT[NTEST] = TEMP2 / max( UNFL, ULP*max( TEMP1, TEMP2 ) );

               } // 680

               NTEST = NTEST + 1;
               slacpy(' ', N, N, V, LDU, A, LDA );
               SRNAMT = 'SSYEVX';
               ssyevx('V', 'I', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK, LWORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSYEVX(V,I,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     RESULT[NTEST+1] = ULPINV;
                     RESULT[NTEST+2] = ULPINV;
                     GO TO 690;
                  }
               }

               // Do tests 31 and 32 (or +54)

               slacpy(' ', N, N, V, LDU, A, LDA );

               ssyt22(1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU, WORK, RESULT( NTEST ) );

               NTEST = NTEST + 2;
               slacpy(' ', N, N, V, LDU, A, LDA );
               SRNAMT = 'SSYEVX';
               ssyevx('N', 'I', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, WORK, LWORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSYEVX(N,I,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     GO TO 690;
                  }
               }

               // Do test 33 (or +54)

               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL );
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL );
               RESULT[NTEST] = ( TEMP1+TEMP2 ) / max( UNFL, ULP*TEMP3 );
               } // 690

               NTEST = NTEST + 1;
               slacpy(' ', N, N, V, LDU, A, LDA );
               SRNAMT = 'SSYEVX';
               ssyevx('V', 'V', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK, LWORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSYEVX(V,V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     RESULT[NTEST+1] = ULPINV;
                     RESULT[NTEST+2] = ULPINV;
                     GO TO 700;
                  }
               }

               // Do tests 34 and 35 (or +54)

               slacpy(' ', N, N, V, LDU, A, LDA );

               ssyt22(1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU, WORK, RESULT( NTEST ) );

               NTEST = NTEST + 2;
               slacpy(' ', N, N, V, LDU, A, LDA );
               SRNAMT = 'SSYEVX';
               ssyevx('N', 'V', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, WORK, LWORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSYEVX(N,V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     GO TO 700;
                  }
               }

               if ( M3 == 0 && N > 0 ) {
                  RESULT[NTEST] = ULPINV;
                  GO TO 700;
               }

               // Do test 36 (or +54)

               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL );
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL );
               if ( N > 0 ) {
                  TEMP3 = max( ( WA1( 1 ) ).abs(), ( WA1( N ) ) ).abs();
               } else {
                  TEMP3 = ZERO;
               }
               RESULT[NTEST] = ( TEMP1+TEMP2 ) / max( UNFL, TEMP3*ULP );

               } // 700

               // 5)      Call SSPEV and SSPEVX.

               slacpy(' ', N, N, V, LDU, A, LDA );

               // Load array WORK with the upper or lower triangular
               // part of the matrix in packed form.

               if ( IUPLO == 1 ) {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 720
                     for (I = 1; I <= J; I++) { // 710
                        WORK[INDX] = A( I, J );
                        INDX = INDX + 1;
                     } // 710
                  } // 720
               } else {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 740
                     for (I = J; I <= N; I++) { // 730
                        WORK[INDX] = A( I, J );
                        INDX = INDX + 1;
                     } // 730
                  } // 740
               }

               NTEST = NTEST + 1;
               SRNAMT = 'SSPEV';
               sspev('V', UPLO, N, WORK, D1, Z, LDU, V, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSPEV(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     RESULT[NTEST+1] = ULPINV;
                     RESULT[NTEST+2] = ULPINV;
                     GO TO 800;
                  }
               }

               // Do tests 37 and 38 (or +54)

               ssyt21(1, UPLO, N, 0, A, LDA, D1, D2, Z, LDU, V, LDU, TAU, WORK, RESULT( NTEST ) );

               if ( IUPLO == 1 ) {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 760
                     for (I = 1; I <= J; I++) { // 750
                        WORK[INDX] = A( I, J );
                        INDX = INDX + 1;
                     } // 750
                  } // 760
               } else {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 780
                     for (I = J; I <= N; I++) { // 770
                        WORK[INDX] = A( I, J );
                        INDX = INDX + 1;
                     } // 770
                  } // 780
               }

               NTEST = NTEST + 2;
               SRNAMT = 'SSPEV';
               sspev('N', UPLO, N, WORK, D3, Z, LDU, V, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSPEV(N,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     GO TO 800;
                  }
               }

               // Do test 39 (or +54)

               TEMP1 = ZERO;
               TEMP2 = ZERO;
               for (J = 1; J <= N; J++) { // 790
                  TEMP1 = max( TEMP1, ( D1( J ) ).abs(), ( D3( J ) ) ).abs();
                  TEMP2 = max( TEMP2, ABS( D1( J )-D3( J ) ) );
               } // 790
               RESULT[NTEST] = TEMP2 / max( UNFL, ULP*max( TEMP1, TEMP2 ) );

               // Load array WORK with the upper or lower triangular part
               // of the matrix in packed form.

               } // 800
               if ( IUPLO == 1 ) {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 820
                     for (I = 1; I <= J; I++) { // 810
                        WORK[INDX] = A( I, J );
                        INDX = INDX + 1;
                     } // 810
                  } // 820
               } else {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 840
                     for (I = J; I <= N; I++) { // 830
                        WORK[INDX] = A( I, J );
                        INDX = INDX + 1;
                     } // 830
                  } // 840
               }

               NTEST = NTEST + 1;

               if ( N > 0 ) {
                  TEMP3 = max( ( D1( 1 ) ).abs(), ( D1( N ) ) ).abs();
                  if ( IL != 1 ) {
                     VL = D1( IL ) - max( HALF*( D1( IL )-D1( IL-1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL );
                  } else if ( N > 0 ) {
                     VL = D1( 1 ) - max( HALF*( D1( N )-D1( 1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL );
                  }
                  if ( IU != N ) {
                     VU = D1( IU ) + max( HALF*( D1( IU+1 )-D1( IU ) ), TEN*ULP*TEMP3, TEN*RTUNFL );
                  } else if ( N > 0 ) {
                     VU = D1( N ) + max( HALF*( D1( N )-D1( 1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL );
                  }
               } else {
                  TEMP3 = ZERO;
                  VL = ZERO;
                  VU = ONE;
               }

               SRNAMT = 'SSPEVX';
               sspevx('V', 'A', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M, WA1, Z, LDU, V, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSPEVX(V,A,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     RESULT[NTEST+1] = ULPINV;
                     RESULT[NTEST+2] = ULPINV;
                     GO TO 900;
                  }
               }

               // Do tests 40 and 41 (or +54)

               ssyt21(1, UPLO, N, 0, A, LDU, WA1, D2, Z, LDU, V, LDU, TAU, WORK, RESULT( NTEST ) );

               NTEST = NTEST + 2;

               if ( IUPLO == 1 ) {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 860
                     for (I = 1; I <= J; I++) { // 850
                        WORK[INDX] = A( I, J );
                        INDX = INDX + 1;
                     } // 850
                  } // 860
               } else {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 880
                     for (I = J; I <= N; I++) { // 870
                        WORK[INDX] = A( I, J );
                        INDX = INDX + 1;
                     } // 870
                  } // 880
               }

               SRNAMT = 'SSPEVX';
               sspevx('N', 'A', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, V, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSPEVX(N,A,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     GO TO 900;
                  }
               }

               // Do test 42 (or +54)

               TEMP1 = ZERO;
               TEMP2 = ZERO;
               for (J = 1; J <= N; J++) { // 890
                  TEMP1 = max( TEMP1, ( WA1( J ) ).abs(), ( WA2( J ) ) ).abs();
                  TEMP2 = max( TEMP2, ABS( WA1( J )-WA2( J ) ) );
               } // 890
               RESULT[NTEST] = TEMP2 / max( UNFL, ULP*max( TEMP1, TEMP2 ) );

               } // 900
               if ( IUPLO == 1 ) {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 920
                     for (I = 1; I <= J; I++) { // 910
                        WORK[INDX] = A( I, J );
                        INDX = INDX + 1;
                     } // 910
                  } // 920
               } else {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 940
                     for (I = J; I <= N; I++) { // 930
                        WORK[INDX] = A( I, J );
                        INDX = INDX + 1;
                     } // 930
                  } // 940
               }

               NTEST = NTEST + 1;

               SRNAMT = 'SSPEVX';
               sspevx('V', 'I', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, V, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSPEVX(V,I,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     RESULT[NTEST+1] = ULPINV;
                     RESULT[NTEST+2] = ULPINV;
                     GO TO 990;
                  }
               }

               // Do tests 43 and 44 (or +54)

               ssyt22(1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU, WORK, RESULT( NTEST ) );

               NTEST = NTEST + 2;

               if ( IUPLO == 1 ) {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 960
                     for (I = 1; I <= J; I++) { // 950
                        WORK[INDX] = A( I, J );
                        INDX = INDX + 1;
                     } // 950
                  } // 960
               } else {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 980
                     for (I = J; I <= N; I++) { // 970
                        WORK[INDX] = A( I, J );
                        INDX = INDX + 1;
                     } // 970
                  } // 980
               }

               SRNAMT = 'SSPEVX';
               sspevx('N', 'I', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, V, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSPEVX(N,I,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     GO TO 990;
                  }
               }

               if ( M3 == 0 && N > 0 ) {
                  RESULT[NTEST] = ULPINV;
                  GO TO 990;
               }

               // Do test 45 (or +54)

               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL );
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL );
               if ( N > 0 ) {
                  TEMP3 = max( ( WA1( 1 ) ).abs(), ( WA1( N ) ) ).abs();
               } else {
                  TEMP3 = ZERO;
               }
               RESULT[NTEST] = ( TEMP1+TEMP2 ) / max( UNFL, TEMP3*ULP );

               } // 990
               if ( IUPLO == 1 ) {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 1010
                     for (I = 1; I <= J; I++) { // 1000
                        WORK[INDX] = A( I, J );
                        INDX = INDX + 1;
                     } // 1000
                  } // 1010
               } else {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 1030
                     for (I = J; I <= N; I++) { // 1020
                        WORK[INDX] = A( I, J );
                        INDX = INDX + 1;
                     } // 1020
                  } // 1030
               }

               NTEST = NTEST + 1;

               SRNAMT = 'SSPEVX';
               sspevx('V', 'V', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, V, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSPEVX(V,V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     RESULT[NTEST+1] = ULPINV;
                     RESULT[NTEST+2] = ULPINV;
                     GO TO 1080;
                  }
               }

               // Do tests 46 and 47 (or +54)

               ssyt22(1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU, WORK, RESULT( NTEST ) );

               NTEST = NTEST + 2;

               if ( IUPLO == 1 ) {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 1050
                     for (I = 1; I <= J; I++) { // 1040
                        WORK[INDX] = A( I, J );
                        INDX = INDX + 1;
                     } // 1040
                  } // 1050
               } else {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 1070
                     for (I = J; I <= N; I++) { // 1060
                        WORK[INDX] = A( I, J );
                        INDX = INDX + 1;
                     } // 1060
                  } // 1070
               }

               SRNAMT = 'SSPEVX';
               sspevx('N', 'V', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, V, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSPEVX(N,V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     GO TO 1080;
                  }
               }

               if ( M3 == 0 && N > 0 ) {
                  RESULT[NTEST] = ULPINV;
                  GO TO 1080;
               }

               // Do test 48 (or +54)

               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL );
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL );
               if ( N > 0 ) {
                  TEMP3 = max( ( WA1( 1 ) ).abs(), ( WA1( N ) ) ).abs();
               } else {
                  TEMP3 = ZERO;
               }
               RESULT[NTEST] = ( TEMP1+TEMP2 ) / max( UNFL, TEMP3*ULP );

               } // 1080

               // 6)      Call SSBEV and SSBEVX.

               if ( JTYPE <= 7 ) {
                  KD = 1;
               } else if ( JTYPE >= 8 && JTYPE <= 15 ) {
                  KD = max( N-1, 0 );
               } else {
                  KD = IHBW;
               }

               // Load array V with the upper or lower triangular part
               // of the matrix in band form.

               if ( IUPLO == 1 ) {
                  for (J = 1; J <= N; J++) { // 1100
                     for (I = max( 1, J-KD ); I <= J; I++) { // 1090
                        V[KD+1+I-J, J] = A( I, J );
                     } // 1090
                  } // 1100
               } else {
                  for (J = 1; J <= N; J++) { // 1120
                     for (I = J; I <= min( N, J+KD ); I++) { // 1110
                        V[1+I-J, J] = A( I, J );
                     } // 1110
                  } // 1120
               }

               NTEST = NTEST + 1;
               SRNAMT = 'SSBEV';
               ssbev('V', UPLO, N, KD, V, LDU, D1, Z, LDU, WORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSBEV(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     RESULT[NTEST+1] = ULPINV;
                     RESULT[NTEST+2] = ULPINV;
                     GO TO 1180;
                  }
               }

               // Do tests 49 and 50 (or ... )

               ssyt21(1, UPLO, N, 0, A, LDA, D1, D2, Z, LDU, V, LDU, TAU, WORK, RESULT( NTEST ) );

               if ( IUPLO == 1 ) {
                  for (J = 1; J <= N; J++) { // 1140
                     for (I = max( 1, J-KD ); I <= J; I++) { // 1130
                        V[KD+1+I-J, J] = A( I, J );
                     } // 1130
                  } // 1140
               } else {
                  for (J = 1; J <= N; J++) { // 1160
                     for (I = J; I <= min( N, J+KD ); I++) { // 1150
                        V[1+I-J, J] = A( I, J );
                     } // 1150
                  } // 1160
               }

               NTEST = NTEST + 2;
               SRNAMT = 'SSBEV';
               ssbev('N', UPLO, N, KD, V, LDU, D3, Z, LDU, WORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSBEV(N,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     GO TO 1180;
                  }
               }

               // Do test 51 (or +54)

               TEMP1 = ZERO;
               TEMP2 = ZERO;
               for (J = 1; J <= N; J++) { // 1170
                  TEMP1 = max( TEMP1, ( D1( J ) ).abs(), ( D3( J ) ) ).abs();
                  TEMP2 = max( TEMP2, ABS( D1( J )-D3( J ) ) );
               } // 1170
               RESULT[NTEST] = TEMP2 / max( UNFL, ULP*max( TEMP1, TEMP2 ) );

               // Load array V with the upper or lower triangular part
               // of the matrix in band form.

               } // 1180
               if ( IUPLO == 1 ) {
                  for (J = 1; J <= N; J++) { // 1200
                     for (I = max( 1, J-KD ); I <= J; I++) { // 1190
                        V[KD+1+I-J, J] = A( I, J );
                     } // 1190
                  } // 1200
               } else {
                  for (J = 1; J <= N; J++) { // 1220
                     for (I = J; I <= min( N, J+KD ); I++) { // 1210
                        V[1+I-J, J] = A( I, J );
                     } // 1210
                  } // 1220
               }

               NTEST = NTEST + 1;
               SRNAMT = 'SSBEVX';
               ssbevx('V', 'A', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL, M, WA2, Z, LDU, WORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSBEVX(V,A,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     RESULT[NTEST+1] = ULPINV;
                     RESULT[NTEST+2] = ULPINV;
                     GO TO 1280;
                  }
               }

               // Do tests 52 and 53 (or +54)

               ssyt21(1, UPLO, N, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU, WORK, RESULT( NTEST ) );

               NTEST = NTEST + 2;

               if ( IUPLO == 1 ) {
                  for (J = 1; J <= N; J++) { // 1240
                     for (I = max( 1, J-KD ); I <= J; I++) { // 1230
                        V[KD+1+I-J, J] = A( I, J );
                     } // 1230
                  } // 1240
               } else {
                  for (J = 1; J <= N; J++) { // 1260
                     for (I = J; I <= min( N, J+KD ); I++) { // 1250
                        V[1+I-J, J] = A( I, J );
                     } // 1250
                  } // 1260
               }

               SRNAMT = 'SSBEVX';
               ssbevx('N', 'A', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, WORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSBEVX(N,A,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     GO TO 1280;
                  }
               }

               // Do test 54 (or +54)

               TEMP1 = ZERO;
               TEMP2 = ZERO;
               for (J = 1; J <= N; J++) { // 1270
                  TEMP1 = max( TEMP1, ( WA2( J ) ).abs(), ( WA3( J ) ) ).abs();
                  TEMP2 = max( TEMP2, ABS( WA2( J )-WA3( J ) ) );
               } // 1270
               RESULT[NTEST] = TEMP2 / max( UNFL, ULP*max( TEMP1, TEMP2 ) );

               } // 1280
               NTEST = NTEST + 1;
               if ( IUPLO == 1 ) {
                  for (J = 1; J <= N; J++) { // 1300
                     for (I = max( 1, J-KD ); I <= J; I++) { // 1290
                        V[KD+1+I-J, J] = A( I, J );
                     } // 1290
                  } // 1300
               } else {
                  for (J = 1; J <= N; J++) { // 1320
                     for (I = J; I <= min( N, J+KD ); I++) { // 1310
                        V[1+I-J, J] = A( I, J );
                     } // 1310
                  } // 1320
               }

               SRNAMT = 'SSBEVX';
               ssbevx('V', 'I', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSBEVX(V,I,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     RESULT[NTEST+1] = ULPINV;
                     RESULT[NTEST+2] = ULPINV;
                     GO TO 1370;
                  }
               }

               // Do tests 55 and 56 (or +54)

               ssyt22(1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU, WORK, RESULT( NTEST ) );

               NTEST = NTEST + 2;

               if ( IUPLO == 1 ) {
                  for (J = 1; J <= N; J++) { // 1340
                     for (I = max( 1, J-KD ); I <= J; I++) { // 1330
                        V[KD+1+I-J, J] = A( I, J );
                     } // 1330
                  } // 1340
               } else {
                  for (J = 1; J <= N; J++) { // 1360
                     for (I = J; I <= min( N, J+KD ); I++) { // 1350
                        V[1+I-J, J] = A( I, J );
                     } // 1350
                  } // 1360
               }

               SRNAMT = 'SSBEVX';
               ssbevx('N', 'I', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, WORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSBEVX(N,I,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     GO TO 1370;
                  }
               }

               // Do test 57 (or +54)

               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL );
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL );
               if ( N > 0 ) {
                  TEMP3 = max( ( WA1( 1 ) ).abs(), ( WA1( N ) ) ).abs();
               } else {
                  TEMP3 = ZERO;
               }
               RESULT[NTEST] = ( TEMP1+TEMP2 ) / max( UNFL, TEMP3*ULP );

               } // 1370
               NTEST = NTEST + 1;
               if ( IUPLO == 1 ) {
                  for (J = 1; J <= N; J++) { // 1390
                     for (I = max( 1, J-KD ); I <= J; I++) { // 1380
                        V[KD+1+I-J, J] = A( I, J );
                     } // 1380
                  } // 1390
               } else {
                  for (J = 1; J <= N; J++) { // 1410
                     for (I = J; I <= min( N, J+KD ); I++) { // 1400
                        V[1+I-J, J] = A( I, J );
                     } // 1400
                  } // 1410
               }

               SRNAMT = 'SSBEVX';
               ssbevx('V', 'V', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSBEVX(V,V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     RESULT[NTEST+1] = ULPINV;
                     RESULT[NTEST+2] = ULPINV;
                     GO TO 1460;
                  }
               }

               // Do tests 58 and 59 (or +54)

               ssyt22(1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU, WORK, RESULT( NTEST ) );

               NTEST = NTEST + 2;

               if ( IUPLO == 1 ) {
                  for (J = 1; J <= N; J++) { // 1430
                     for (I = max( 1, J-KD ); I <= J; I++) { // 1420
                        V[KD+1+I-J, J] = A( I, J );
                     } // 1420
                  } // 1430
               } else {
                  for (J = 1; J <= N; J++) { // 1450
                     for (I = J; I <= min( N, J+KD ); I++) { // 1440
                        V[1+I-J, J] = A( I, J );
                     } // 1440
                  } // 1450
               }

               SRNAMT = 'SSBEVX';
               ssbevx('N', 'V', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, WORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSBEVX(N,V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     GO TO 1460;
                  }
               }

               if ( M3 == 0 && N > 0 ) {
                  RESULT[NTEST] = ULPINV;
                  GO TO 1460;
               }

               // Do test 60 (or +54)

               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL );
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL );
               if ( N > 0 ) {
                  TEMP3 = max( ( WA1( 1 ) ).abs(), ( WA1( N ) ) ).abs();
               } else {
                  TEMP3 = ZERO;
               }
               RESULT[NTEST] = ( TEMP1+TEMP2 ) / max( UNFL, TEMP3*ULP );

               } // 1460

               // 7)      Call SSYEVD

               slacpy(' ', N, N, A, LDA, V, LDU );

               NTEST = NTEST + 1;
               SRNAMT = 'SSYEVD';
               ssyevd('V', UPLO, N, A, LDU, D1, WORK, LWEDC, IWORK, LIWEDC, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSYEVD(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     RESULT[NTEST+1] = ULPINV;
                     RESULT[NTEST+2] = ULPINV;
                     GO TO 1480;
                  }
               }

               // Do tests 61 and 62 (or +54)

               ssyt21(1, UPLO, N, 0, V, LDU, D1, D2, A, LDU, Z, LDU, TAU, WORK, RESULT( NTEST ) );

               slacpy(' ', N, N, V, LDU, A, LDA );

               NTEST = NTEST + 2;
               SRNAMT = 'SSYEVD';
               ssyevd('N', UPLO, N, A, LDU, D3, WORK, LWEDC, IWORK, LIWEDC, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSYEVD(N,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     GO TO 1480;
                  }
               }

               // Do test 63 (or +54)

               TEMP1 = ZERO;
               TEMP2 = ZERO;
               for (J = 1; J <= N; J++) { // 1470
                  TEMP1 = max( TEMP1, ( D1( J ) ).abs(), ( D3( J ) ) ).abs();
                  TEMP2 = max( TEMP2, ABS( D1( J )-D3( J ) ) );
               } // 1470
               RESULT[NTEST] = TEMP2 / max( UNFL, ULP*max( TEMP1, TEMP2 ) );

               } // 1480

               // 8)      Call SSPEVD.

               slacpy(' ', N, N, V, LDU, A, LDA );

               // Load array WORK with the upper or lower triangular
               // part of the matrix in packed form.

               if ( IUPLO == 1 ) {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 1500
                     for (I = 1; I <= J; I++) { // 1490
                        WORK[INDX] = A( I, J );
                        INDX = INDX + 1;
                     } // 1490
                  } // 1500
               } else {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 1520
                     for (I = J; I <= N; I++) { // 1510
                        WORK[INDX] = A( I, J );
                        INDX = INDX + 1;
                     } // 1510
                  } // 1520
               }

               NTEST = NTEST + 1;
               SRNAMT = 'SSPEVD';
               sspevd('V', UPLO, N, WORK, D1, Z, LDU, WORK( INDX ), LWEDC-INDX+1, IWORK, LIWEDC, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSPEVD(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     RESULT[NTEST+1] = ULPINV;
                     RESULT[NTEST+2] = ULPINV;
                     GO TO 1580;
                  }
               }

               // Do tests 64 and 65 (or +54)

               ssyt21(1, UPLO, N, 0, A, LDA, D1, D2, Z, LDU, V, LDU, TAU, WORK, RESULT( NTEST ) );

               if ( IUPLO == 1 ) {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 1540
                     for (I = 1; I <= J; I++) { // 1530

                        WORK[INDX] = A( I, J );
                        INDX = INDX + 1;
                     } // 1530
                  } // 1540
               } else {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 1560
                     for (I = J; I <= N; I++) { // 1550
                        WORK[INDX] = A( I, J );
                        INDX = INDX + 1;
                     } // 1550
                  } // 1560
               }

               NTEST = NTEST + 2;
               SRNAMT = 'SSPEVD';
               sspevd('N', UPLO, N, WORK, D3, Z, LDU, WORK( INDX ), LWEDC-INDX+1, IWORK, LIWEDC, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSPEVD(N,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     GO TO 1580;
                  }
               }

               // Do test 66 (or +54)

               TEMP1 = ZERO;
               TEMP2 = ZERO;
               for (J = 1; J <= N; J++) { // 1570
                  TEMP1 = max( TEMP1, ( D1( J ) ).abs(), ( D3( J ) ) ).abs();
                  TEMP2 = max( TEMP2, ABS( D1( J )-D3( J ) ) );
               } // 1570
               RESULT[NTEST] = TEMP2 / max( UNFL, ULP*max( TEMP1, TEMP2 ) );
               } // 1580

               // 9)      Call SSBEVD.

               if ( JTYPE <= 7 ) {
                  KD = 1;
               } else if ( JTYPE >= 8 && JTYPE <= 15 ) {
                  KD = max( N-1, 0 );
               } else {
                  KD = IHBW;
               }

               // Load array V with the upper or lower triangular part
               // of the matrix in band form.

               if ( IUPLO == 1 ) {
                  for (J = 1; J <= N; J++) { // 1600
                     for (I = max( 1, J-KD ); I <= J; I++) { // 1590
                        V[KD+1+I-J, J] = A( I, J );
                     } // 1590
                  } // 1600
               } else {
                  for (J = 1; J <= N; J++) { // 1620
                     for (I = J; I <= min( N, J+KD ); I++) { // 1610
                        V[1+I-J, J] = A( I, J );
                     } // 1610
                  } // 1620
               }

               NTEST = NTEST + 1;
               SRNAMT = 'SSBEVD';
               ssbevd('V', UPLO, N, KD, V, LDU, D1, Z, LDU, WORK, LWEDC, IWORK, LIWEDC, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSBEVD(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     RESULT[NTEST+1] = ULPINV;
                     RESULT[NTEST+2] = ULPINV;
                     GO TO 1680;
                  }
               }

               // Do tests 67 and 68 (or +54)

               ssyt21(1, UPLO, N, 0, A, LDA, D1, D2, Z, LDU, V, LDU, TAU, WORK, RESULT( NTEST ) );

               if ( IUPLO == 1 ) {
                  for (J = 1; J <= N; J++) { // 1640
                     for (I = max( 1, J-KD ); I <= J; I++) { // 1630
                        V[KD+1+I-J, J] = A( I, J );
                     } // 1630
                  } // 1640
               } else {
                  for (J = 1; J <= N; J++) { // 1660
                     for (I = J; I <= min( N, J+KD ); I++) { // 1650
                        V[1+I-J, J] = A( I, J );
                     } // 1650
                  } // 1660
               }

               NTEST = NTEST + 2;
               SRNAMT = 'SSBEVD';
               ssbevd('N', UPLO, N, KD, V, LDU, D3, Z, LDU, WORK, LWEDC, IWORK, LIWEDC, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSBEVD(N,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     GO TO 1680;
                  }
               }

               // Do test 69 (or +54)

               TEMP1 = ZERO;
               TEMP2 = ZERO;
               for (J = 1; J <= N; J++) { // 1670
                  TEMP1 = max( TEMP1, ( D1( J ) ).abs(), ( D3( J ) ) ).abs();
                  TEMP2 = max( TEMP2, ABS( D1( J )-D3( J ) ) );
               } // 1670
               RESULT[NTEST] = TEMP2 / max( UNFL, ULP*max( TEMP1, TEMP2 ) );

               } // 1680


               slacpy(' ', N, N, A, LDA, V, LDU );
               NTEST = NTEST + 1;
               SRNAMT = 'SSYEVR';
               ssyevr('V', 'A', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M, WA1, Z, LDU, IWORK, WORK, LWORK, IWORK(2*N+1), LIWORK-2*N, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSYEVR(V,A,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     RESULT[NTEST+1] = ULPINV;
                     RESULT[NTEST+2] = ULPINV;
                     GO TO 1700;
                  }
               }

               // Do tests 70 and 71 (or ... )

               slacpy(' ', N, N, V, LDU, A, LDA );

               ssyt21(1, UPLO, N, 0, A, LDU, WA1, D2, Z, LDU, V, LDU, TAU, WORK, RESULT( NTEST ) );

               NTEST = NTEST + 2;
               SRNAMT = 'SSYEVR';
               ssyevr('N', 'A', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, IWORK, WORK, LWORK, IWORK(2*N+1), LIWORK-2*N, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSYEVR(N,A,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     GO TO 1700;
                  }
               }

               // Do test 72 (or ... )

               TEMP1 = ZERO;
               TEMP2 = ZERO;
               for (J = 1; J <= N; J++) { // 1690
                  TEMP1 = max( TEMP1, ( WA1( J ) ).abs(), ( WA2( J ) ) ).abs();
                  TEMP2 = max( TEMP2, ABS( WA1( J )-WA2( J ) ) );
               } // 1690
               RESULT[NTEST] = TEMP2 / max( UNFL, ULP*max( TEMP1, TEMP2 ) );

               } // 1700

               NTEST = NTEST + 1;
               slacpy(' ', N, N, V, LDU, A, LDA );
               SRNAMT = 'SSYEVR';
               ssyevr('V', 'I', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, IWORK, WORK, LWORK, IWORK(2*N+1), LIWORK-2*N, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSYEVR(V,I,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     RESULT[NTEST+1] = ULPINV;
                     RESULT[NTEST+2] = ULPINV;
                     GO TO 1710;
                  }
               }

               // Do tests 73 and 74 (or +54)

               slacpy(' ', N, N, V, LDU, A, LDA );

               ssyt22(1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU, WORK, RESULT( NTEST ) );

               NTEST = NTEST + 2;
               slacpy(' ', N, N, V, LDU, A, LDA );
               SRNAMT = 'SSYEVR';
               ssyevr('N', 'I', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, IWORK, WORK, LWORK, IWORK(2*N+1), LIWORK-2*N, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSYEVR(N,I,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     GO TO 1710;
                  }
               }

               // Do test 75 (or +54)

               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL );
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL );
               RESULT[NTEST] = ( TEMP1+TEMP2 ) / max( UNFL, ULP*TEMP3 );
               } // 1710

               NTEST = NTEST + 1;
               slacpy(' ', N, N, V, LDU, A, LDA );
               SRNAMT = 'SSYEVR';
               ssyevr('V', 'V', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, IWORK, WORK, LWORK, IWORK(2*N+1), LIWORK-2*N, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSYEVR(V,V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     RESULT[NTEST+1] = ULPINV;
                     RESULT[NTEST+2] = ULPINV;
                     GO TO 700;
                  }
               }

               // Do tests 76 and 77 (or +54)

               slacpy(' ', N, N, V, LDU, A, LDA );

               ssyt22(1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU, WORK, RESULT( NTEST ) );

               NTEST = NTEST + 2;
               slacpy(' ', N, N, V, LDU, A, LDA );
               SRNAMT = 'SSYEVR';
               ssyevr('N', 'V', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, IWORK, WORK, LWORK, IWORK(2*N+1), LIWORK-2*N, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'SSYEVR(N,V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  if ( IINFO < 0 ) {
                     return;
                  } else {
                     RESULT[NTEST] = ULPINV;
                     GO TO 700;
                  }
               }

               if ( M3 == 0 && N > 0 ) {
                  RESULT[NTEST] = ULPINV;
                  GO TO 700;
               }

               // Do test 78 (or +54)

               TEMP1 = SSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL );
               TEMP2 = SSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL );
               if ( N > 0 ) {
                  TEMP3 = max( ( WA1( 1 ) ).abs(), ( WA1( N ) ) ).abs();
               } else {
                  TEMP3 = ZERO;
               }
               RESULT[NTEST] = ( TEMP1+TEMP2 ) / max( UNFL, TEMP3*ULP );

               slacpy(' ', N, N, V, LDU, A, LDA );

            } // 1720

            // End of Loop -- Check for RESULT(j) > THRESH

            NTESTT = NTESTT + NTEST;

            slafts('SST', N, N, JTYPE, NTEST, RESULT, IOLDSD, THRESH, NOUNIT, NERRS );

         } // 1730
      } // 1740

      // Summary

      alasvm('SST', NOUNIT, NERRS, NTESTT, 0 );

 9999 FORMAT( ' SDRVST: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' );

      return;
      }