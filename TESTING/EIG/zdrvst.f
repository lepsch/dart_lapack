      SUBROUTINE ZDRVST( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH, NOUNIT, A, LDA, D1, D2, D3, WA1, WA2, WA3, U, LDU, V, TAU, Z, WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, RESULT, INFO );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDU, LIWORK, LRWORK, LWORK, NOUNIT, NSIZES, NTYPES;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                ISEED( 4 ), IWORK( * ), NN( * );
      double             D1( * ), D2( * ), D3( * ), RESULT( * ), RWORK( * ), WA1( * ), WA2( * ), WA3( * );
      COMPLEX*16         A( LDA, * ), TAU( * ), U( LDU, * ), V( LDU, * ), WORK( * ), Z( LDU, * );
      // ..

*  =====================================================================


      // .. Parameters ..
      double             ZERO, ONE, TWO, TEN;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0, TEN = 10.0 ;
      double             HALF;
      const              HALF = ONE / TWO ;
      COMPLEX*16         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      int                MAXTYP;
      const              MAXTYP = 18 ;
      // ..
      // .. Local Scalars ..
      bool               BADNN;
      String             UPLO;
      int                I, IDIAG, IHBW, IINFO, IL, IMODE, INDWRK, INDX, IROW, ITEMP, ITYPE, IU, IUPLO, J, J1, J2, JCOL, JSIZE, JTYPE, KD, LGN, LIWEDC, LRWEDC, LWEDC, M, M2, M3, MTYPES, N, NERRS, NMATS, NMAX, NTEST, NTESTT;
      double             ABSTOL, ANINV, ANORM, COND, OVFL, RTOVFL, RTUNFL, TEMP1, TEMP2, TEMP3, ULP, ULPINV, UNFL, VL, VU;
      // ..
      // .. Local Arrays ..
      int                IDUMMA( 1 ), IOLDSD( 4 ), ISEED2( 4 ), ISEED3( 4 ), KMAGN( MAXTYP ), KMODE( MAXTYP ), KTYPE( MAXTYP );
      // ..
      // .. External Functions ..
      double             DLAMCH, DLARND, DSXT1;
      // EXTERNAL DLAMCH, DLARND, DSXT1
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALASVM, DLAFTS, XERBLA, ZHBEV, ZHBEVD, ZHBEVX, ZHEEV, ZHEEVD, ZHEEVR, ZHEEVX, ZHET21, ZHET22, ZHPEV, ZHPEVD, ZHPEVX, ZLACPY, ZLASET, ZLATMR, ZLATMS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, INT, LOG, MAX, MIN, SQRT
      // ..
      // .. Data statements ..
      DATA               KTYPE / 1, 2, 5*4, 5*5, 3*8, 3*9 /;
      DATA               KMAGN / 2*1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3, 1, 2, 3 /;
      DATA               KMODE / 2*0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0, 4, 4, 4 /;
      // ..
      // .. Executable Statements ..

      // 1)      Check for errors

      NTESTT = 0;
      INFO = 0;

      BADNN = false;
      NMAX = 1;
      for (J = 1; J <= NSIZES; J++) { // 10
         NMAX = MAX( NMAX, NN( J ) );
         IF( NN( J ) < 0 ) BADNN = true;
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
      } else if ( 2*MAX( 2, NMAX )**2 > LWORK ) {
         INFO = -22;
      }

      if ( INFO != 0 ) {
         xerbla('ZDRVST', -INFO );
         RETURN;
      }

      // Quick return if nothing to do

      if (NSIZES == 0 || NTYPES == 0) RETURN;

      // More Important constants

      UNFL = DLAMCH( 'Safe minimum' );
      OVFL = DLAMCH( 'Overflow' );
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' );
      ULPINV = ONE / ULP;
      RTUNFL = SQRT( UNFL );
      RTOVFL = SQRT( OVFL );

      // Loop over sizes, types

      for (I = 1; I <= 4; I++) { // 20
         ISEED2( I ) = ISEED( I );
         ISEED3( I ) = ISEED( I );
      } // 20

      NERRS = 0;
      NMATS = 0;

      for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) { // 1220
         N = NN( JSIZE );
         if ( N > 0 ) {
            LGN = INT( LOG( DBLE( N ) ) / LOG( TWO ) );
            if (2**LGN < N) LGN = LGN + 1             IF( 2**LGN < N ) LGN = LGN + 1;
            LWEDC = MAX( 2*N+N*N, 2*N*N );
            LRWEDC = 1 + 4*N + 2*N*LGN + 3*N**2;
            LIWEDC = 3 + 5*N;
         } else {
            LWEDC = 2;
            LRWEDC = 8;
            LIWEDC = 8;
         }
         ANINV = ONE / DBLE( MAX( 1, N ) );

         if ( NSIZES != 1 ) {
            MTYPES = MIN( MAXTYP, NTYPES );
         } else {
            MTYPES = MIN( MAXTYP+1, NTYPES );
         }

         for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) { // 1210
            IF( !DOTYPE( JTYPE ) ) GO TO 1210;
            NMATS = NMATS + 1;
            NTEST = 0;

            for (J = 1; J <= 4; J++) { // 30
               IOLDSD( J ) = ISEED( J );
            } // 30

            // 2)      Compute "A"

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
            // =9                      band Hermitian, w/ eigenvalues

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

            zlaset('Full', LDA, N, CZERO, CZERO, A, LDA );
            IINFO = 0;
            COND = ULPINV;

            // Special Matrices -- Identity & Jordan block

                    // Zero

            if ( ITYPE == 1 ) {
               IINFO = 0;

            } else if ( ITYPE == 2 ) {

               // Identity

               for (JCOL = 1; JCOL <= N; JCOL++) { // 80
                  A( JCOL, JCOL ) = ANORM;
               } // 80

            } else if ( ITYPE == 4 ) {

               // Diagonal Matrix, [Eigen]values Specified

               zlatms(N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, 0, 0, 'N', A, LDA, WORK, IINFO );

            } else if ( ITYPE == 5 ) {

               // Hermitian, eigenvalues specified

               zlatms(N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, N, N, 'N', A, LDA, WORK, IINFO );

            } else if ( ITYPE == 7 ) {

               // Diagonal, random eigenvalues

               zlatmr(N, N, 'S', ISEED, 'H', WORK, 6, ONE, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, 0, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else if ( ITYPE == 8 ) {

               // Hermitian, random eigenvalues

               zlatmr(N, N, 'S', ISEED, 'H', WORK, 6, ONE, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else if ( ITYPE == 9 ) {

               // Hermitian banded, eigenvalues specified

               IHBW = INT( ( N-1 )*DLARND( 1, ISEED3 ) );
               zlatms(N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, IHBW, IHBW, 'Z', U, LDU, WORK, IINFO );

               // Store as dense matrix for most routines.

               zlaset('Full', LDA, N, CZERO, CZERO, A, LDA );
               for (IDIAG = -IHBW; IDIAG <= IHBW; IDIAG++) { // 100
                  IROW = IHBW - IDIAG + 1;
                  J1 = MAX( 1, IDIAG+1 );
                  J2 = MIN( N, N+IDIAG );
                  for (J = J1; J <= J2; J++) { // 90
                     I = J - IDIAG;
                     A( I, J ) = U( IROW, J );
                  } // 90
               } // 100
            } else {
               IINFO = 1;
            }

            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'Generator', IINFO, N, JTYPE, IOLDSD;
               INFO = ABS( IINFO );
               RETURN;
            }

            } // 110

            ABSTOL = UNFL + UNFL;
            if ( N <= 1 ) {
               IL = 1;
               IU = N;
            } else {
               IL = 1 + INT( ( N-1 )*DLARND( 1, ISEED2 ) );
               IU = 1 + INT( ( N-1 )*DLARND( 1, ISEED2 ) );
               if ( IL > IU ) {
                  ITEMP = IL;
                  IL = IU;
                  IU = ITEMP;
               }
            }

            // Perform tests storing upper or lower triangular
            // part of matrix.

            for (IUPLO = 0; IUPLO <= 1; IUPLO++) { // 1200
               if ( IUPLO == 0 ) {
                  UPLO = 'L';
               } else {
                  UPLO = 'U';
               }

               // Call ZHEEVD and CHEEVX.

               zlacpy(' ', N, N, A, LDA, V, LDU );

               NTEST = NTEST + 1;
               zheevd('V', UPLO, N, A, LDU, D1, WORK, LWEDC, RWORK, LRWEDC, IWORK, LIWEDC, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'ZHEEVD(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     RESULT( NTEST+1 ) = ULPINV;
                     RESULT( NTEST+2 ) = ULPINV;
                     GO TO 130;
                  }
               }

               // Do tests 1 and 2.

               zhet21(1, UPLO, N, 0, V, LDU, D1, D2, A, LDU, Z, LDU, TAU, WORK, RWORK, RESULT( NTEST ) );

               zlacpy(' ', N, N, V, LDU, A, LDA );

               NTEST = NTEST + 2;
               zheevd('N', UPLO, N, A, LDU, D3, WORK, LWEDC, RWORK, LRWEDC, IWORK, LIWEDC, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'ZHEEVD(N,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     GO TO 130;
                  }
               }

               // Do test 3.

               TEMP1 = ZERO;
               TEMP2 = ZERO;
               for (J = 1; J <= N; J++) { // 120
                  TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D3( J ) ) );
                  TEMP2 = MAX( TEMP2, ABS( D1( J )-D3( J ) ) );
               } // 120
               RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) );

               } // 130
               zlacpy(' ', N, N, V, LDU, A, LDA );

               NTEST = NTEST + 1;

               if ( N > 0 ) {
                  TEMP3 = MAX( ABS( D1( 1 ) ), ABS( D1( N ) ) );
                  if ( IL != 1 ) {
                     VL = D1( IL ) - MAX( HALF*( D1( IL )-D1( IL-1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL );
                  } else if ( N > 0 ) {
                     VL = D1( 1 ) - MAX( HALF*( D1( N )-D1( 1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL );
                  }
                  if ( IU != N ) {
                     VU = D1( IU ) + MAX( HALF*( D1( IU+1 )-D1( IU ) ), TEN*ULP*TEMP3, TEN*RTUNFL );
                  } else if ( N > 0 ) {
                     VU = D1( N ) + MAX( HALF*( D1( N )-D1( 1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL );
                  }
               } else {
                  TEMP3 = ZERO;
                  VL = ZERO;
                  VU = ONE;
               }

               zheevx('V', 'A', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M, WA1, Z, LDU, WORK, LWORK, RWORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'ZHEEVX(V,A,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     RESULT( NTEST+1 ) = ULPINV;
                     RESULT( NTEST+2 ) = ULPINV;
                     GO TO 150;
                  }
               }

               // Do tests 4 and 5.

               zlacpy(' ', N, N, V, LDU, A, LDA );

               zhet21(1, UPLO, N, 0, A, LDU, WA1, D2, Z, LDU, V, LDU, TAU, WORK, RWORK, RESULT( NTEST ) );

               NTEST = NTEST + 2;
               zheevx('N', 'A', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK, LWORK, RWORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'ZHEEVX(N,A,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     GO TO 150;
                  }
               }

               // Do test 6.

               TEMP1 = ZERO;
               TEMP2 = ZERO;
               for (J = 1; J <= N; J++) { // 140
                  TEMP1 = MAX( TEMP1, ABS( WA1( J ) ), ABS( WA2( J ) ) );
                  TEMP2 = MAX( TEMP2, ABS( WA1( J )-WA2( J ) ) );
               } // 140
               RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) );

               } // 150
               zlacpy(' ', N, N, V, LDU, A, LDA );

               NTEST = NTEST + 1;

               zheevx('V', 'I', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK, LWORK, RWORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'ZHEEVX(V,I,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     GO TO 160;
                  }
               }

               // Do tests 7 and 8.

               zlacpy(' ', N, N, V, LDU, A, LDA );

               zhet22(1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU, WORK, RWORK, RESULT( NTEST ) );

               NTEST = NTEST + 2;

               zheevx('N', 'I', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, WORK, LWORK, RWORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'ZHEEVX(N,I,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     GO TO 160;
                  }
               }

               // Do test 9.

               TEMP1 = DSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL );
               TEMP2 = DSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL );
               if ( N > 0 ) {
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) );
               } else {
                  TEMP3 = ZERO;
               }
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) / MAX( UNFL, TEMP3*ULP );

               } // 160
               zlacpy(' ', N, N, V, LDU, A, LDA );

               NTEST = NTEST + 1;

               zheevx('V', 'V', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK, LWORK, RWORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'ZHEEVX(V,V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     GO TO 170;
                  }
               }

               // Do tests 10 and 11.

               zlacpy(' ', N, N, V, LDU, A, LDA );

               zhet22(1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU, WORK, RWORK, RESULT( NTEST ) );

               NTEST = NTEST + 2;

               zheevx('N', 'V', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, WORK, LWORK, RWORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'ZHEEVX(N,V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     GO TO 170;
                  }
               }

               if ( M3 == 0 && N > 0 ) {
                  RESULT( NTEST ) = ULPINV;
                  GO TO 170;
               }

               // Do test 12.

               TEMP1 = DSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL );
               TEMP2 = DSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL );
               if ( N > 0 ) {
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) );
               } else {
                  TEMP3 = ZERO;
               }
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) / MAX( UNFL, TEMP3*ULP );

               } // 170

               // Call ZHPEVD and CHPEVX.

               zlacpy(' ', N, N, V, LDU, A, LDA );

               // Load array WORK with the upper or lower triangular
               // part of the matrix in packed form.

               if ( IUPLO == 1 ) {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 190
                     for (I = 1; I <= J; I++) { // 180
                        WORK( INDX ) = A( I, J );
                        INDX = INDX + 1;
                     } // 180
                  } // 190
               } else {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 210
                     for (I = J; I <= N; I++) { // 200
                        WORK( INDX ) = A( I, J );
                        INDX = INDX + 1;
                     } // 200
                  } // 210
               }

               NTEST = NTEST + 1;
               INDWRK = N*( N+1 ) / 2 + 1;
               zhpevd('V', UPLO, N, WORK, D1, Z, LDU, WORK( INDWRK ), LWEDC, RWORK, LRWEDC, IWORK, LIWEDC, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'ZHPEVD(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     RESULT( NTEST+1 ) = ULPINV;
                     RESULT( NTEST+2 ) = ULPINV;
                     GO TO 270;
                  }
               }

               // Do tests 13 and 14.

               zhet21(1, UPLO, N, 0, A, LDA, D1, D2, Z, LDU, V, LDU, TAU, WORK, RWORK, RESULT( NTEST ) );

               if ( IUPLO == 1 ) {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 230
                     for (I = 1; I <= J; I++) { // 220
                        WORK( INDX ) = A( I, J );
                        INDX = INDX + 1;
                     } // 220
                  } // 230
               } else {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 250
                     for (I = J; I <= N; I++) { // 240
                        WORK( INDX ) = A( I, J );
                        INDX = INDX + 1;
                     } // 240
                  } // 250
               }

               NTEST = NTEST + 2;
               INDWRK = N*( N+1 ) / 2 + 1;
               zhpevd('N', UPLO, N, WORK, D3, Z, LDU, WORK( INDWRK ), LWEDC, RWORK, LRWEDC, IWORK, LIWEDC, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'ZHPEVD(N,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     GO TO 270;
                  }
               }

               // Do test 15.

               TEMP1 = ZERO;
               TEMP2 = ZERO;
               for (J = 1; J <= N; J++) { // 260
                  TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D3( J ) ) );
                  TEMP2 = MAX( TEMP2, ABS( D1( J )-D3( J ) ) );
               } // 260
               RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) );

               // Load array WORK with the upper or lower triangular part
               // of the matrix in packed form.

               } // 270
               if ( IUPLO == 1 ) {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 290
                     for (I = 1; I <= J; I++) { // 280
                        WORK( INDX ) = A( I, J );
                        INDX = INDX + 1;
                     } // 280
                  } // 290
               } else {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 310
                     for (I = J; I <= N; I++) { // 300
                        WORK( INDX ) = A( I, J );
                        INDX = INDX + 1;
                     } // 300
                  } // 310
               }

               NTEST = NTEST + 1;

               if ( N > 0 ) {
                  TEMP3 = MAX( ABS( D1( 1 ) ), ABS( D1( N ) ) );
                  if ( IL != 1 ) {
                     VL = D1( IL ) - MAX( HALF*( D1( IL )-D1( IL-1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL );
                  } else if ( N > 0 ) {
                     VL = D1( 1 ) - MAX( HALF*( D1( N )-D1( 1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL );
                  }
                  if ( IU != N ) {
                     VU = D1( IU ) + MAX( HALF*( D1( IU+1 )-D1( IU ) ), TEN*ULP*TEMP3, TEN*RTUNFL );
                  } else if ( N > 0 ) {
                     VU = D1( N ) + MAX( HALF*( D1( N )-D1( 1 ) ), TEN*ULP*TEMP3, TEN*RTUNFL );
                  }
               } else {
                  TEMP3 = ZERO;
                  VL = ZERO;
                  VU = ONE;
               }

               zhpevx('V', 'A', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M, WA1, Z, LDU, V, RWORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'ZHPEVX(V,A,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     RESULT( NTEST+1 ) = ULPINV;
                     RESULT( NTEST+2 ) = ULPINV;
                     GO TO 370;
                  }
               }

               // Do tests 16 and 17.

               zhet21(1, UPLO, N, 0, A, LDU, WA1, D2, Z, LDU, V, LDU, TAU, WORK, RWORK, RESULT( NTEST ) );

               NTEST = NTEST + 2;

               if ( IUPLO == 1 ) {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 330
                     for (I = 1; I <= J; I++) { // 320
                        WORK( INDX ) = A( I, J );
                        INDX = INDX + 1;
                     } // 320
                  } // 330
               } else {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 350
                     for (I = J; I <= N; I++) { // 340
                        WORK( INDX ) = A( I, J );
                        INDX = INDX + 1;
                     } // 340
                  } // 350
               }

               zhpevx('N', 'A', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, V, RWORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'ZHPEVX(N,A,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     GO TO 370;
                  }
               }

               // Do test 18.

               TEMP1 = ZERO;
               TEMP2 = ZERO;
               for (J = 1; J <= N; J++) { // 360
                  TEMP1 = MAX( TEMP1, ABS( WA1( J ) ), ABS( WA2( J ) ) );
                  TEMP2 = MAX( TEMP2, ABS( WA1( J )-WA2( J ) ) );
               } // 360
               RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) );

               } // 370
               NTEST = NTEST + 1;
               if ( IUPLO == 1 ) {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 390
                     for (I = 1; I <= J; I++) { // 380
                        WORK( INDX ) = A( I, J );
                        INDX = INDX + 1;
                     } // 380
                  } // 390
               } else {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 410
                     for (I = J; I <= N; I++) { // 400
                        WORK( INDX ) = A( I, J );
                        INDX = INDX + 1;
                     } // 400
                  } // 410
               }

               zhpevx('V', 'I', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, V, RWORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'ZHPEVX(V,I,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     RESULT( NTEST+1 ) = ULPINV;
                     RESULT( NTEST+2 ) = ULPINV;
                     GO TO 460;
                  }
               }

               // Do tests 19 and 20.

               zhet22(1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU, WORK, RWORK, RESULT( NTEST ) );

               NTEST = NTEST + 2;

               if ( IUPLO == 1 ) {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 430
                     for (I = 1; I <= J; I++) { // 420
                        WORK( INDX ) = A( I, J );
                        INDX = INDX + 1;
                     } // 420
                  } // 430
               } else {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 450
                     for (I = J; I <= N; I++) { // 440
                        WORK( INDX ) = A( I, J );
                        INDX = INDX + 1;
                     } // 440
                  } // 450
               }

               zhpevx('N', 'I', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, V, RWORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'ZHPEVX(N,I,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     GO TO 460;
                  }
               }

               // Do test 21.

               TEMP1 = DSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL );
               TEMP2 = DSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL );
               if ( N > 0 ) {
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) );
               } else {
                  TEMP3 = ZERO;
               }
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) / MAX( UNFL, TEMP3*ULP );

               } // 460
               NTEST = NTEST + 1;
               if ( IUPLO == 1 ) {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 480
                     for (I = 1; I <= J; I++) { // 470
                        WORK( INDX ) = A( I, J );
                        INDX = INDX + 1;
                     } // 470
                  } // 480
               } else {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 500
                     for (I = J; I <= N; I++) { // 490
                        WORK( INDX ) = A( I, J );
                        INDX = INDX + 1;
                     } // 490
                  } // 500
               }

               zhpevx('V', 'V', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, V, RWORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'ZHPEVX(V,V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     RESULT( NTEST+1 ) = ULPINV;
                     RESULT( NTEST+2 ) = ULPINV;
                     GO TO 550;
                  }
               }

               // Do tests 22 and 23.

               zhet22(1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU, WORK, RWORK, RESULT( NTEST ) );

               NTEST = NTEST + 2;

               if ( IUPLO == 1 ) {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 520
                     for (I = 1; I <= J; I++) { // 510
                        WORK( INDX ) = A( I, J );
                        INDX = INDX + 1;
                     } // 510
                  } // 520
               } else {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 540
                     for (I = J; I <= N; I++) { // 530
                        WORK( INDX ) = A( I, J );
                        INDX = INDX + 1;
                     } // 530
                  } // 540
               }

               zhpevx('N', 'V', UPLO, N, WORK, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, V, RWORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'ZHPEVX(N,V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     GO TO 550;
                  }
               }

               if ( M3 == 0 && N > 0 ) {
                  RESULT( NTEST ) = ULPINV;
                  GO TO 550;
               }

               // Do test 24.

               TEMP1 = DSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL );
               TEMP2 = DSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL );
               if ( N > 0 ) {
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) );
               } else {
                  TEMP3 = ZERO;
               }
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) / MAX( UNFL, TEMP3*ULP );

               } // 550

               // Call ZHBEVD and CHBEVX.

               if ( JTYPE <= 7 ) {
                  KD = 0;
               } else if ( JTYPE >= 8 && JTYPE <= 15 ) {
                  KD = MAX( N-1, 0 );
               } else {
                  KD = IHBW;
               }

               // Load array V with the upper or lower triangular part
               // of the matrix in band form.

               if ( IUPLO == 1 ) {
                  for (J = 1; J <= N; J++) { // 570
                     DO 560 I = MAX( 1, J-KD ), J;
                        V( KD+1+I-J, J ) = A( I, J );
                     } // 560
                  } // 570
               } else {
                  for (J = 1; J <= N; J++) { // 590
                     DO 580 I = J, MIN( N, J+KD );
                        V( 1+I-J, J ) = A( I, J );
                     } // 580
                  } // 590
               }

               NTEST = NTEST + 1;
               zhbevd('V', UPLO, N, KD, V, LDU, D1, Z, LDU, WORK, LWEDC, RWORK, LRWEDC, IWORK, LIWEDC, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9998 )'ZHBEVD(V,' // UPLO // ')', IINFO, N, KD, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     RESULT( NTEST+1 ) = ULPINV;
                     RESULT( NTEST+2 ) = ULPINV;
                     GO TO 650;
                  }
               }

               // Do tests 25 and 26.

               zhet21(1, UPLO, N, 0, A, LDA, D1, D2, Z, LDU, V, LDU, TAU, WORK, RWORK, RESULT( NTEST ) );

               if ( IUPLO == 1 ) {
                  for (J = 1; J <= N; J++) { // 610
                     DO 600 I = MAX( 1, J-KD ), J;
                        V( KD+1+I-J, J ) = A( I, J );
                     } // 600
                  } // 610
               } else {
                  for (J = 1; J <= N; J++) { // 630
                     DO 620 I = J, MIN( N, J+KD );
                        V( 1+I-J, J ) = A( I, J );
                     } // 620
                  } // 630
               }

               NTEST = NTEST + 2;
               zhbevd('N', UPLO, N, KD, V, LDU, D3, Z, LDU, WORK, LWEDC, RWORK, LRWEDC, IWORK, LIWEDC, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9998 )'ZHBEVD(N,' // UPLO // ')', IINFO, N, KD, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     GO TO 650;
                  }
               }

               // Do test 27.

               TEMP1 = ZERO;
               TEMP2 = ZERO;
               for (J = 1; J <= N; J++) { // 640
                  TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D3( J ) ) );
                  TEMP2 = MAX( TEMP2, ABS( D1( J )-D3( J ) ) );
               } // 640
               RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) );

               // Load array V with the upper or lower triangular part
               // of the matrix in band form.

               } // 650
               if ( IUPLO == 1 ) {
                  for (J = 1; J <= N; J++) { // 670
                     DO 660 I = MAX( 1, J-KD ), J;
                        V( KD+1+I-J, J ) = A( I, J );
                     } // 660
                  } // 670
               } else {
                  for (J = 1; J <= N; J++) { // 690
                     DO 680 I = J, MIN( N, J+KD );
                        V( 1+I-J, J ) = A( I, J );
                     } // 680
                  } // 690
               }

               NTEST = NTEST + 1;
               zhbevx('V', 'A', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL, M, WA1, Z, LDU, WORK, RWORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'ZHBEVX(V,A,' // UPLO // ')', IINFO, N, KD, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     RESULT( NTEST+1 ) = ULPINV;
                     RESULT( NTEST+2 ) = ULPINV;
                     GO TO 750;
                  }
               }

               // Do tests 28 and 29.

               zhet21(1, UPLO, N, 0, A, LDU, WA1, D2, Z, LDU, V, LDU, TAU, WORK, RWORK, RESULT( NTEST ) );

               NTEST = NTEST + 2;

               if ( IUPLO == 1 ) {
                  for (J = 1; J <= N; J++) { // 710
                     DO 700 I = MAX( 1, J-KD ), J;
                        V( KD+1+I-J, J ) = A( I, J );
                     } // 700
                  } // 710
               } else {
                  for (J = 1; J <= N; J++) { // 730
                     DO 720 I = J, MIN( N, J+KD );
                        V( 1+I-J, J ) = A( I, J );
                     } // 720
                  } // 730
               }

               zhbevx('N', 'A', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK, RWORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9998 )'ZHBEVX(N,A,' // UPLO // ')', IINFO, N, KD, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     GO TO 750;
                  }
               }

               // Do test 30.

               TEMP1 = ZERO;
               TEMP2 = ZERO;
               for (J = 1; J <= N; J++) { // 740
                  TEMP1 = MAX( TEMP1, ABS( WA1( J ) ), ABS( WA2( J ) ) );
                  TEMP2 = MAX( TEMP2, ABS( WA1( J )-WA2( J ) ) );
               } // 740
               RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) );

               // Load array V with the upper or lower triangular part
               // of the matrix in band form.

               } // 750
               NTEST = NTEST + 1;
               if ( IUPLO == 1 ) {
                  for (J = 1; J <= N; J++) { // 770
                     DO 760 I = MAX( 1, J-KD ), J;
                        V( KD+1+I-J, J ) = A( I, J );
                     } // 760
                  } // 770
               } else {
                  for (J = 1; J <= N; J++) { // 790
                     DO 780 I = J, MIN( N, J+KD );
                        V( 1+I-J, J ) = A( I, J );
                     } // 780
                  } // 790
               }

               zhbevx('V', 'I', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK, RWORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9998 )'ZHBEVX(V,I,' // UPLO // ')', IINFO, N, KD, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     RESULT( NTEST+1 ) = ULPINV;
                     RESULT( NTEST+2 ) = ULPINV;
                     GO TO 840;
                  }
               }

               // Do tests 31 and 32.

               zhet22(1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU, WORK, RWORK, RESULT( NTEST ) );

               NTEST = NTEST + 2;

               if ( IUPLO == 1 ) {
                  for (J = 1; J <= N; J++) { // 810
                     DO 800 I = MAX( 1, J-KD ), J;
                        V( KD+1+I-J, J ) = A( I, J );
                     } // 800
                  } // 810
               } else {
                  for (J = 1; J <= N; J++) { // 830
                     DO 820 I = J, MIN( N, J+KD );
                        V( 1+I-J, J ) = A( I, J );
                     } // 820
                  } // 830
               }
               zhbevx('N', 'I', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, WORK, RWORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9998 )'ZHBEVX(N,I,' // UPLO // ')', IINFO, N, KD, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     GO TO 840;
                  }
               }

               // Do test 33.

               TEMP1 = DSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL );
               TEMP2 = DSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL );
               if ( N > 0 ) {
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) );
               } else {
                  TEMP3 = ZERO;
               }
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) / MAX( UNFL, TEMP3*ULP );

               // Load array V with the upper or lower triangular part
               // of the matrix in band form.

               } // 840
               NTEST = NTEST + 1;
               if ( IUPLO == 1 ) {
                  for (J = 1; J <= N; J++) { // 860
                     DO 850 I = MAX( 1, J-KD ), J;
                        V( KD+1+I-J, J ) = A( I, J );
                     } // 850
                  } // 860
               } else {
                  for (J = 1; J <= N; J++) { // 880
                     DO 870 I = J, MIN( N, J+KD );
                        V( 1+I-J, J ) = A( I, J );
                     } // 870
                  } // 880
               }
               zhbevx('V', 'V', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, WORK, RWORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9998 )'ZHBEVX(V,V,' // UPLO // ')', IINFO, N, KD, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     RESULT( NTEST+1 ) = ULPINV;
                     RESULT( NTEST+2 ) = ULPINV;
                     GO TO 930;
                  }
               }

               // Do tests 34 and 35.

               zhet22(1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU, WORK, RWORK, RESULT( NTEST ) );

               NTEST = NTEST + 2;

               if ( IUPLO == 1 ) {
                  for (J = 1; J <= N; J++) { // 900
                     DO 890 I = MAX( 1, J-KD ), J;
                        V( KD+1+I-J, J ) = A( I, J );
                     } // 890
                  } // 900
               } else {
                  for (J = 1; J <= N; J++) { // 920
                     DO 910 I = J, MIN( N, J+KD );
                        V( 1+I-J, J ) = A( I, J );
                     } // 910
                  } // 920
               }
               zhbevx('N', 'V', UPLO, N, KD, V, LDU, U, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, WORK, RWORK, IWORK, IWORK( 5*N+1 ), IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9998 )'ZHBEVX(N,V,' // UPLO // ')', IINFO, N, KD, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     GO TO 930;
                  }
               }

               if ( M3 == 0 && N > 0 ) {
                  RESULT( NTEST ) = ULPINV;
                  GO TO 930;
               }

               // Do test 36.

               TEMP1 = DSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL );
               TEMP2 = DSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL );
               if ( N > 0 ) {
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) );
               } else {
                  TEMP3 = ZERO;
               }
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) / MAX( UNFL, TEMP3*ULP );

               } // 930

               // Call ZHEEV

               zlacpy(' ', N, N, A, LDA, V, LDU );

               NTEST = NTEST + 1;
               zheev('V', UPLO, N, A, LDU, D1, WORK, LWORK, RWORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'ZHEEV(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     RESULT( NTEST+1 ) = ULPINV;
                     RESULT( NTEST+2 ) = ULPINV;
                     GO TO 950;
                  }
               }

               // Do tests 37 and 38

               zhet21(1, UPLO, N, 0, V, LDU, D1, D2, A, LDU, Z, LDU, TAU, WORK, RWORK, RESULT( NTEST ) );

               zlacpy(' ', N, N, V, LDU, A, LDA );

               NTEST = NTEST + 2;
               zheev('N', UPLO, N, A, LDU, D3, WORK, LWORK, RWORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'ZHEEV(N,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     GO TO 950;
                  }
               }

               // Do test 39

               TEMP1 = ZERO;
               TEMP2 = ZERO;
               for (J = 1; J <= N; J++) { // 940
                  TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D3( J ) ) );
                  TEMP2 = MAX( TEMP2, ABS( D1( J )-D3( J ) ) );
               } // 940
               RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) );

               } // 950

               zlacpy(' ', N, N, V, LDU, A, LDA );

               // Call ZHPEV

               // Load array WORK with the upper or lower triangular
               // part of the matrix in packed form.

               if ( IUPLO == 1 ) {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 970
                     for (I = 1; I <= J; I++) { // 960
                        WORK( INDX ) = A( I, J );
                        INDX = INDX + 1;
                     } // 960
                  } // 970
               } else {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 990
                     for (I = J; I <= N; I++) { // 980
                        WORK( INDX ) = A( I, J );
                        INDX = INDX + 1;
                     } // 980
                  } // 990
               }

               NTEST = NTEST + 1;
               INDWRK = N*( N+1 ) / 2 + 1;
               zhpev('V', UPLO, N, WORK, D1, Z, LDU, WORK( INDWRK ), RWORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'ZHPEV(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     RESULT( NTEST+1 ) = ULPINV;
                     RESULT( NTEST+2 ) = ULPINV;
                     GO TO 1050;
                  }
               }

               // Do tests 40 and 41.

               zhet21(1, UPLO, N, 0, A, LDA, D1, D2, Z, LDU, V, LDU, TAU, WORK, RWORK, RESULT( NTEST ) );

               if ( IUPLO == 1 ) {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 1010
                     for (I = 1; I <= J; I++) { // 1000
                        WORK( INDX ) = A( I, J );
                        INDX = INDX + 1;
                     } // 1000
                  } // 1010
               } else {
                  INDX = 1;
                  for (J = 1; J <= N; J++) { // 1030
                     for (I = J; I <= N; I++) { // 1020
                        WORK( INDX ) = A( I, J );
                        INDX = INDX + 1;
                     } // 1020
                  } // 1030
               }

               NTEST = NTEST + 2;
               INDWRK = N*( N+1 ) / 2 + 1;
               zhpev('N', UPLO, N, WORK, D3, Z, LDU, WORK( INDWRK ), RWORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'ZHPEV(N,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     GO TO 1050;
                  }
               }

               // Do test 42

               TEMP1 = ZERO;
               TEMP2 = ZERO;
               for (J = 1; J <= N; J++) { // 1040
                  TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D3( J ) ) );
                  TEMP2 = MAX( TEMP2, ABS( D1( J )-D3( J ) ) );
               } // 1040
               RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) );

               } // 1050

               // Call ZHBEV

               if ( JTYPE <= 7 ) {
                  KD = 0;
               } else if ( JTYPE >= 8 && JTYPE <= 15 ) {
                  KD = MAX( N-1, 0 );
               } else {
                  KD = IHBW;
               }

               // Load array V with the upper or lower triangular part
               // of the matrix in band form.

               if ( IUPLO == 1 ) {
                  for (J = 1; J <= N; J++) { // 1070
                     DO 1060 I = MAX( 1, J-KD ), J;
                        V( KD+1+I-J, J ) = A( I, J );
                     } // 1060
                  } // 1070
               } else {
                  for (J = 1; J <= N; J++) { // 1090
                     DO 1080 I = J, MIN( N, J+KD );
                        V( 1+I-J, J ) = A( I, J );
                     } // 1080
                  } // 1090
               }

               NTEST = NTEST + 1;
               zhbev('V', UPLO, N, KD, V, LDU, D1, Z, LDU, WORK, RWORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9998 )'ZHBEV(V,' // UPLO // ')', IINFO, N, KD, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     RESULT( NTEST+1 ) = ULPINV;
                     RESULT( NTEST+2 ) = ULPINV;
                     GO TO 1140;
                  }
               }

               // Do tests 43 and 44.

               zhet21(1, UPLO, N, 0, A, LDA, D1, D2, Z, LDU, V, LDU, TAU, WORK, RWORK, RESULT( NTEST ) );

               if ( IUPLO == 1 ) {
                  for (J = 1; J <= N; J++) { // 1110
                     DO 1100 I = MAX( 1, J-KD ), J;
                        V( KD+1+I-J, J ) = A( I, J );
                     } // 1100
                  } // 1110
               } else {
                  for (J = 1; J <= N; J++) { // 1130
                     DO 1120 I = J, MIN( N, J+KD );
                        V( 1+I-J, J ) = A( I, J );
                     } // 1120
                  } // 1130
               }

               NTEST = NTEST + 2;
               zhbev('N', UPLO, N, KD, V, LDU, D3, Z, LDU, WORK, RWORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9998 )'ZHBEV(N,' // UPLO // ')', IINFO, N, KD, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     GO TO 1140;
                  }
               }

               } // 1140

               // Do test 45.

               TEMP1 = ZERO;
               TEMP2 = ZERO;
               for (J = 1; J <= N; J++) { // 1150
                  TEMP1 = MAX( TEMP1, ABS( D1( J ) ), ABS( D3( J ) ) );
                  TEMP2 = MAX( TEMP2, ABS( D1( J )-D3( J ) ) );
               } // 1150
               RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) );

               zlacpy(' ', N, N, A, LDA, V, LDU );
               NTEST = NTEST + 1;
               zheevr('V', 'A', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M, WA1, Z, LDU, IWORK, WORK, LWORK, RWORK, LRWORK, IWORK( 2*N+1 ), LIWORK-2*N, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'ZHEEVR(V,A,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     RESULT( NTEST+1 ) = ULPINV;
                     RESULT( NTEST+2 ) = ULPINV;
                     GO TO 1170;
                  }
               }

               // Do tests 45 and 46 (or ... )

               zlacpy(' ', N, N, V, LDU, A, LDA );

               zhet21(1, UPLO, N, 0, A, LDU, WA1, D2, Z, LDU, V, LDU, TAU, WORK, RWORK, RESULT( NTEST ) );

               NTEST = NTEST + 2;
               zheevr('N', 'A', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, IWORK, WORK, LWORK, RWORK, LRWORK, IWORK( 2*N+1 ), LIWORK-2*N, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'ZHEEVR(N,A,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     GO TO 1170;
                  }
               }

               // Do test 47 (or ... )

               TEMP1 = ZERO;
               TEMP2 = ZERO;
               for (J = 1; J <= N; J++) { // 1160
                  TEMP1 = MAX( TEMP1, ABS( WA1( J ) ), ABS( WA2( J ) ) );
                  TEMP2 = MAX( TEMP2, ABS( WA1( J )-WA2( J ) ) );
               } // 1160
               RESULT( NTEST ) = TEMP2 / MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) );

               } // 1170

               NTEST = NTEST + 1;
               zlacpy(' ', N, N, V, LDU, A, LDA );
               zheevr('V', 'I', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, IWORK, WORK, LWORK, RWORK, LRWORK, IWORK( 2*N+1 ), LIWORK-2*N, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'ZHEEVR(V,I,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     RESULT( NTEST+1 ) = ULPINV;
                     RESULT( NTEST+2 ) = ULPINV;
                     GO TO 1180;
                  }
               }

               // Do tests 48 and 49 (or +??)

               zlacpy(' ', N, N, V, LDU, A, LDA );

               zhet22(1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU, WORK, RWORK, RESULT( NTEST ) );

               NTEST = NTEST + 2;
               zlacpy(' ', N, N, V, LDU, A, LDA );
               zheevr('N', 'I', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, IWORK, WORK, LWORK, RWORK, LRWORK, IWORK( 2*N+1 ), LIWORK-2*N, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'ZHEEVR(N,I,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     GO TO 1180;
                  }
               }

               // Do test 50 (or +??)

               TEMP1 = DSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL );
               TEMP2 = DSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL );
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) / MAX( UNFL, ULP*TEMP3 );
               } // 1180

               NTEST = NTEST + 1;
               zlacpy(' ', N, N, V, LDU, A, LDA );
               zheevr('V', 'V', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M2, WA2, Z, LDU, IWORK, WORK, LWORK, RWORK, LRWORK, IWORK( 2*N+1 ), LIWORK-2*N, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'ZHEEVR(V,V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     RESULT( NTEST+1 ) = ULPINV;
                     RESULT( NTEST+2 ) = ULPINV;
                     GO TO 1190;
                  }
               }

               // Do tests 51 and 52 (or +??)

               zlacpy(' ', N, N, V, LDU, A, LDA );

               zhet22(1, UPLO, N, M2, 0, A, LDU, WA2, D2, Z, LDU, V, LDU, TAU, WORK, RWORK, RESULT( NTEST ) );

               NTEST = NTEST + 2;
               zlacpy(' ', N, N, V, LDU, A, LDA );
               zheevr('N', 'V', UPLO, N, A, LDU, VL, VU, IL, IU, ABSTOL, M3, WA3, Z, LDU, IWORK, WORK, LWORK, RWORK, LRWORK, IWORK( 2*N+1 ), LIWORK-2*N, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9999 )'ZHEEVR(N,V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                  INFO = ABS( IINFO );
                  if ( IINFO < 0 ) {
                     RETURN;
                  } else {
                     RESULT( NTEST ) = ULPINV;
                     GO TO 1190;
                  }
               }

               if ( M3 == 0 && N > 0 ) {
                  RESULT( NTEST ) = ULPINV;
                  GO TO 1190;
               }

               // Do test 52 (or +??)

               TEMP1 = DSXT1( 1, WA2, M2, WA3, M3, ABSTOL, ULP, UNFL );
               TEMP2 = DSXT1( 1, WA3, M3, WA2, M2, ABSTOL, ULP, UNFL );
               if ( N > 0 ) {
                  TEMP3 = MAX( ABS( WA1( 1 ) ), ABS( WA1( N ) ) );
               } else {
                  TEMP3 = ZERO;
               }
               RESULT( NTEST ) = ( TEMP1+TEMP2 ) / MAX( UNFL, TEMP3*ULP );

               zlacpy(' ', N, N, V, LDU, A, LDA );




               // Load array V with the upper or lower triangular part
               // of the matrix in band form.

               } // 1190

            } // 1200

            // End of Loop -- Check for RESULT(j) > THRESH

            NTESTT = NTESTT + NTEST;
            dlafts('ZST', N, N, JTYPE, NTEST, RESULT, IOLDSD, THRESH, NOUNIT, NERRS );

         } // 1210
      } // 1220

      // Summary

      alasvm('ZST', NOUNIT, NERRS, NTESTT, 0 );

 9999 FORMAT( ' ZDRVST: ', A, ' returned INFO=', I6, / 9X, 'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' );
 9998 FORMAT( ' ZDRVST: ', A, ' returned INFO=', I6, / 9X, 'N=', I6, ', KD=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' );

      RETURN;

      // End of ZDRVST

      }
