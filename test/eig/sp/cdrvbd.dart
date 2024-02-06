      void cdrvbd(NSIZES, MM, NN, NTYPES, DOTYPE, ISEED, THRESH, A, LDA, U, LDU, VT, LDVT, ASAV, USAV, VTSAV, S, SSAV, E, WORK, LWORK, RWORK, IWORK, NOUNIT, INFO ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // IMPLICIT NONE
      int                INFO, LDA, LDU, LDVT, LWORK, NOUNIT, NSIZES, NTYPES;
      double               THRESH;
      bool               DOTYPE( * );
      int                ISEED( 4 ), IWORK( * ), MM( * ), NN( * );
      double               E( * ), RWORK( * ), S( * ), SSAV( * );
      Complex            A( LDA, * ), ASAV( LDA, * ), U( LDU, * ), USAV( LDU, * ), VT( LDVT, * ), VTSAV( LDVT, * ), WORK( * );
      // ..

      double               ZERO, ONE, TWO, HALF;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0, HALF = 0.5 ;
      Complex            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      int                MAXTYP;
      const              MAXTYP = 5 ;
      bool               BADMM, BADNN;
      String             JOBQ, JOBU, JOBVT, RANGE;
      int                I, IINFO, IJQ, IJU, IJVT, IL, IU, ITEMP, IWSPC, IWTMP, J, JSIZE, JTYPE, LSWORK, M, MINWRK, MMAX, MNMAX, MNMIN, MTYPES, N, NERRS, NFAIL, NMAX, NS, NSI, NSV, NTEST, NTESTF, NTESTT, LRWORK;
      double               ANORM, DIF, DIV, OVFL, RTUNFL, ULP, ULPINV, UNFL, VL, VU;
      // ..
      // .. Local Scalars for CGESVDQ ..
      int                LIWORK, NUMRANK;
      String             CJOB( 4 ), CJOBR( 3 ), CJOBV( 2 );
      int                IOLDSD( 4 ), ISEED2( 4 );
      double               RESULT( 39 );
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH, SLARND;
      // EXTERNAL SLAMCH, SLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALASVM, XERBLA, CBDT01, CBDT05, CGESDD, CGESVD, CGESVDQ, CGESVJ, CGEJSV, CGESVDX, CLACPY, CLASET, CLATMS, CUNT01, CUNT03
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, REAL, MAX, MIN
      // ..
      // .. Scalars in Common ..
      String            srnamc.SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Data statements ..
      const CJOB = [ 'N', 'O', 'S', 'A' ];
      const CJOBR = [ 'A', 'V', 'I' ];
      const CJOBV = [ 'N', 'V' ];

      // Check for errors

      INFO = 0;

      // Important constants

      NERRS = 0;
      NTESTT = 0;
      NTESTF = 0;
      BADMM = false;
      BADNN = false;
      MMAX = 1;
      NMAX = 1;
      MNMAX = 1;
      MINWRK = 1;
      for (J = 1; J <= NSIZES; J++) { // 10
         MMAX = max( MMAX, MM( J ) );
         if( MM( J ) < 0 ) BADMM = true;
         NMAX = max( NMAX, NN( J ) );
         if( NN( J ) < 0 ) BADNN = true;
         MNMAX = max( MNMAX, min( MM( J ), NN( J ) ) );
         MINWRK = max( MINWRK, max( 3*min( MM( J ), NN( J ) )+max( MM( J ), NN( J ) )**2, 5*min( MM( J ), NN( J ) ), 3*max( MM( J ), NN( J ) ) ) );
      } // 10

      // Check for errors

      if ( NSIZES < 0 ) {
         INFO = -1;
      } else if ( BADMM ) {
         INFO = -2;
      } else if ( BADNN ) {
         INFO = -3;
      } else if ( NTYPES < 0 ) {
         INFO = -4;
      } else if ( LDA < max( 1, MMAX ) ) {
         INFO = -10;
      } else if ( LDU < max( 1, MMAX ) ) {
         INFO = -12;
      } else if ( LDVT < max( 1, NMAX ) ) {
         INFO = -14;
      } else if ( MINWRK > LWORK ) {
         INFO = -21;
      }

      if ( INFO != 0 ) {
         xerbla('CDRVBD', -INFO );
         return;
      }

      // Quick return if nothing to do

      if (NSIZES == 0 || NTYPES == 0) return;

      // More Important constants

      UNFL = SLAMCH( 'S' );
      OVFL = ONE / UNFL;
      ULP = SLAMCH( 'E' );
      ULPINV = ONE / ULP;
      RTUNFL = sqrt( UNFL );

      // Loop over sizes, types

      NERRS = 0;

      for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) { // 310
         M = MM( JSIZE );
         N = NN( JSIZE );
         MNMIN = min( M, N );

         if ( NSIZES != 1 ) {
            MTYPES = min( MAXTYP, NTYPES );
         } else {
            MTYPES = min( MAXTYP+1, NTYPES );
         }

         for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) { // 300
            if( !DOTYPE( JTYPE ) ) GO TO 300;
            NTEST = 0;

            for (J = 1; J <= 4; J++) { // 20
               IOLDSD[J] = ISEED( J );
            } // 20

            // Compute "A"

            if (MTYPES > MAXTYP) GO TO 50;

            if ( JTYPE == 1 ) {

               // Zero matrix

               claset('Full', M, N, CZERO, CZERO, A, LDA );
               for (I = 1; I <= min( M, N ); I++) { // 30
                  S[I] = ZERO;
               } // 30

            } else if ( JTYPE == 2 ) {

               // Identity matrix

               claset('Full', M, N, CZERO, CONE, A, LDA );
               for (I = 1; I <= min( M, N ); I++) { // 40
                  S[I] = ONE;
               } // 40

            } else {

               // (Scaled) random matrix

               if (JTYPE == 3) ANORM = ONE;
               if( JTYPE == 4 ) ANORM = UNFL / ULP;
               IF( JTYPE == 5 ) ANORM = OVFL*ULP;
               clatms(M, N, 'U', ISEED, 'N', S, 4, double( MNMIN ), ANORM, M-1, N-1, 'N', A, LDA, WORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9996 )'Generator', IINFO, M, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  return;
               }
            }

            } // 50
            clacpy('F', M, N, A, LDA, ASAV, LDA );

            // Do for minimal and adequate (for blocking) workspace

            for (IWSPC = 1; IWSPC <= 4; IWSPC++) { // 290

               // Test for CGESVD

               IWTMP = 2*min( M, N )+max( M, N );
               LSWORK = IWTMP + ( IWSPC-1 )*( LWORK-IWTMP ) / 3;
               LSWORK = min( LSWORK, LWORK );
               LSWORK = max( LSWORK, 1 );
               if (IWSPC == 4) LSWORK = LWORK;

               for (J = 1; J <= 35; J++) { // 60
                  RESULT[J] = -ONE;
               } // 60

               // Factorize A

               if (IWSPC > 1) clacpy( 'F', M, N, ASAV, LDA, A, LDA );
              srnamc.SRNAMT = 'CGESVD';
               cgesvd('A', 'A', M, N, A, LDA, SSAV, USAV, LDU, VTSAV, LDVT, WORK, LSWORK, RWORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9995 )'GESVD', IINFO, M, N, JTYPE, LSWORK, IOLDSD;
                  INFO = ( IINFO ).abs();
                  return;
               }

               // Do tests 1--4

               cbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RWORK, RESULT( 1 ) );
               if ( M != 0 && N != 0 ) {
                  cunt01('Columns', MNMIN, M, USAV, LDU, WORK, LWORK, RWORK, RESULT( 2 ) );
                  cunt01('Rows', MNMIN, N, VTSAV, LDVT, WORK, LWORK, RWORK, RESULT( 3 ) );
               }
               RESULT[4] = 0;
               for (I = 1; I <= MNMIN - 1; I++) { // 70
                  if[SSAV( I ) < SSAV( I+1 ) ) RESULT( 4] = ULPINV;
                  IF[SSAV( I ) < ZERO ) RESULT( 4] = ULPINV;
               } // 70
               if ( MNMIN >= 1 ) {
                  if[SSAV( MNMIN ) < ZERO ) RESULT( 4] = ULPINV;
               }

               // Do partial SVDs, comparing to SSAV, USAV, and VTSAV

               RESULT[5] = ZERO;
               RESULT[6] = ZERO;
               RESULT[7] = ZERO;
               for (IJU = 0; IJU <= 3; IJU++) { // 100
                  for (IJVT = 0; IJVT <= 3; IJVT++) { // 90
                     if( ( IJU == 3 && IJVT == 3 ) || ( IJU == 1 && IJVT == 1 ) )GO TO 90;
                     JOBU = CJOB( IJU+1 );
                     JOBVT = CJOB( IJVT+1 );
                     clacpy('F', M, N, ASAV, LDA, A, LDA );
                    srnamc.SRNAMT = 'CGESVD';
                     cgesvd(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LSWORK, RWORK, IINFO );

                     // Compare U

                     DIF = ZERO;
                     if ( M > 0 && N > 0 ) {
                        if ( IJU == 1 ) {
                           cunt03('C', M, MNMIN, M, MNMIN, USAV, LDU, A, LDA, WORK, LWORK, RWORK, DIF, IINFO );
                        } else if ( IJU == 2 ) {
                           cunt03('C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, RWORK, DIF, IINFO );
                        } else if ( IJU == 3 ) {
                           cunt03('C', M, M, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, RWORK, DIF, IINFO );
                        }
                     }
                     RESULT[5] = max( RESULT( 5 ), DIF );

                     // Compare VT

                     DIF = ZERO;
                     if ( M > 0 && N > 0 ) {
                        if ( IJVT == 1 ) {
                           cunt03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, A, LDA, WORK, LWORK, RWORK, DIF, IINFO );
                        } else if ( IJVT == 2 ) {
                           cunt03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, RWORK, DIF, IINFO );
                        } else if ( IJVT == 3 ) {
                           cunt03('R', N, N, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, RWORK, DIF, IINFO );
                        }
                     }
                     RESULT[6] = max( RESULT( 6 ), DIF );

                     // Compare S

                     DIF = ZERO;
                     DIV = max( double( MNMIN )*ULP*S( 1 ), SLAMCH( 'Safe minimum' ) );
                     for (I = 1; I <= MNMIN - 1; I++) { // 80
                        if( SSAV( I ) < SSAV( I+1 ) ) DIF = ULPINV;
                        IF( SSAV( I ) < ZERO ) DIF = ULPINV;
                        DIF = max( DIF, ABS( SSAV( I )-S( I ) ) / DIV );
                     } // 80
                     RESULT[7] = max( RESULT( 7 ), DIF );
                  } // 90
               } // 100

               // Test for CGESDD

               IWTMP = 2*MNMIN*MNMIN + 2*MNMIN + max( M, N );
               LSWORK = IWTMP + ( IWSPC-1 )*( LWORK-IWTMP ) / 3;
               LSWORK = min( LSWORK, LWORK );
               LSWORK = max( LSWORK, 1 );
               if (IWSPC == 4) LSWORK = LWORK;

               // Factorize A

               clacpy('F', M, N, ASAV, LDA, A, LDA );
              srnamc.SRNAMT = 'CGESDD';
               cgesdd('A', M, N, A, LDA, SSAV, USAV, LDU, VTSAV, LDVT, WORK, LSWORK, RWORK, IWORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9995 )'GESDD', IINFO, M, N, JTYPE, LSWORK, IOLDSD;
                  INFO = ( IINFO ).abs();
                  return;
               }

               // Do tests 1--4

               cbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RWORK, RESULT( 8 ) );
               if ( M != 0 && N != 0 ) {
                  cunt01('Columns', MNMIN, M, USAV, LDU, WORK, LWORK, RWORK, RESULT( 9 ) );
                  cunt01('Rows', MNMIN, N, VTSAV, LDVT, WORK, LWORK, RWORK, RESULT( 10 ) );
               }
               RESULT[11] = 0;
               for (I = 1; I <= MNMIN - 1; I++) { // 110
                  if[SSAV( I ) < SSAV( I+1 ) ) RESULT( 11] = ULPINV;
                  IF[SSAV( I ) < ZERO ) RESULT( 11] = ULPINV;
               } // 110
               if ( MNMIN >= 1 ) {
                  if[SSAV( MNMIN ) < ZERO ) RESULT( 11] = ULPINV;
               }

               // Do partial SVDs, comparing to SSAV, USAV, and VTSAV

               RESULT[12] = ZERO;
               RESULT[13] = ZERO;
               RESULT[14] = ZERO;
               for (IJQ = 0; IJQ <= 2; IJQ++) { // 130
                  JOBQ = CJOB( IJQ+1 );
                  clacpy('F', M, N, ASAV, LDA, A, LDA );
                 srnamc.SRNAMT = 'CGESDD';
                  cgesdd(JOBQ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LSWORK, RWORK, IWORK, IINFO );

                  // Compare U

                  DIF = ZERO;
                  if ( M > 0 && N > 0 ) {
                     if ( IJQ == 1 ) {
                        if ( M >= N ) {
                           cunt03('C', M, MNMIN, M, MNMIN, USAV, LDU, A, LDA, WORK, LWORK, RWORK, DIF, IINFO );
                        } else {
                           cunt03('C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, RWORK, DIF, IINFO );
                        }
                     } else if ( IJQ == 2 ) {
                        cunt03('C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, RWORK, DIF, IINFO );
                     }
                  }
                  RESULT[12] = max( RESULT( 12 ), DIF );

                  // Compare VT

                  DIF = ZERO;
                  if ( M > 0 && N > 0 ) {
                     if ( IJQ == 1 ) {
                        if ( M >= N ) {
                           cunt03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, RWORK, DIF, IINFO );
                        } else {
                           cunt03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, A, LDA, WORK, LWORK, RWORK, DIF, IINFO );
                        }
                     } else if ( IJQ == 2 ) {
                        cunt03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, RWORK, DIF, IINFO );
                     }
                  }
                  RESULT[13] = max( RESULT( 13 ), DIF );

                  // Compare S

                  DIF = ZERO;
                  DIV = max( double( MNMIN )*ULP*S( 1 ), SLAMCH( 'Safe minimum' ) );
                  for (I = 1; I <= MNMIN - 1; I++) { // 120
                     if( SSAV( I ) < SSAV( I+1 ) ) DIF = ULPINV;
                     IF( SSAV( I ) < ZERO ) DIF = ULPINV;
                     DIF = max( DIF, ABS( SSAV( I )-S( I ) ) / DIV );
                  } // 120
                  RESULT[14] = max( RESULT( 14 ), DIF );
               } // 130


               // Test CGESVDQ
               // Note: CGESVDQ only works for M >= N

               RESULT[36] = ZERO;
               RESULT[37] = ZERO;
               RESULT[38] = ZERO;
               RESULT[39] = ZERO;

               if ( M >= N ) {
                  IWTMP = 2*MNMIN*MNMIN + 2*MNMIN + max( M, N );
                  LSWORK = IWTMP + ( IWSPC-1 )*( LWORK-IWTMP ) / 3;
                  LSWORK = min( LSWORK, LWORK );
                  LSWORK = max( LSWORK, 1 );
                  if (IWSPC == 4) LSWORK = LWORK;

                  clacpy('F', M, N, ASAV, LDA, A, LDA );
                 srnamc.SRNAMT = 'CGESVDQ';

                  LRWORK = max(2, M, 5*N);
                  LIWORK = max( N, 1 );
                  cgesvdq('H', 'N', 'N', 'A', 'A',  M, N, A, LDA, SSAV, USAV, LDU, VTSAV, LDVT, NUMRANK, IWORK, LIWORK, WORK, LWORK, RWORK, LRWORK, IINFO );

                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9995 )'CGESVDQ', IINFO, M, N, JTYPE, LSWORK, IOLDSD;
                     INFO = ( IINFO ).abs();
                     return;
                  }

                  // Do tests 36--39

                  cbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RWORK, RESULT( 36 ) );
                  if ( M != 0 && N != 0 ) {
                     cunt01('Columns', M, M, USAV, LDU, WORK, LWORK, RWORK, RESULT( 37 ) );
                     cunt01('Rows', N, N, VTSAV, LDVT, WORK, LWORK, RWORK, RESULT( 38 ) );
                  }
                  RESULT[39] = ZERO;
                  for (I = 1; I <= MNMIN - 1; I++) { // 199
                     if[SSAV( I ) < SSAV( I+1 ) ) RESULT( 39] = ULPINV;
                     IF[SSAV( I ) < ZERO ) RESULT( 39] = ULPINV;
                  } // 199
                  if ( MNMIN >= 1 ) {
                     if[SSAV( MNMIN ) < ZERO ) RESULT( 39] = ULPINV;
                  }
               }

               // Test CGESVJ
               // Note: CGESVJ only works for M >= N

               RESULT[15] = ZERO;
               RESULT[16] = ZERO;
               RESULT[17] = ZERO;
               RESULT[18] = ZERO;

               if ( M >= N ) {
                  IWTMP = 2*MNMIN*MNMIN + 2*MNMIN + max( M, N );
                  LSWORK = IWTMP + ( IWSPC-1 )*( LWORK-IWTMP ) / 3;
                  LSWORK = min( LSWORK, LWORK );
                  LSWORK = max( LSWORK, 1 );
                  LRWORK = max(6,N);
                  if (IWSPC == 4) LSWORK = LWORK;

                  clacpy('F', M, N, ASAV, LDA, USAV, LDA );
                 srnamc.SRNAMT = 'CGESVJ';
                  cgesvj('G', 'U', 'V', M, N, USAV, LDA, SSAV, 0, A, LDVT, WORK, LWORK, RWORK, LRWORK, IINFO );

                  // CGESVJ returns V not VH

                  for (J = 1; J <= N; J++) {
                     for (I = 1; I <= N; I++) {
                        VTSAV[J][I] = CONJG (A(I,J));
                     }
                  }

                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9995 )'GESVJ', IINFO, M, N, JTYPE, LSWORK, IOLDSD;
                     INFO = ( IINFO ).abs();
                     return;
                  }

                  // Do tests 15--18

                  cbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RWORK, RESULT( 15 ) );
                  if ( M != 0 && N != 0 ) {
                     cunt01('Columns', M, M, USAV, LDU, WORK, LWORK, RWORK, RESULT( 16 ) );
                     cunt01('Rows', N, N, VTSAV, LDVT, WORK, LWORK, RWORK, RESULT( 17 ) );
                  }
                  RESULT[18] = ZERO;
                  for (I = 1; I <= MNMIN - 1; I++) { // 131
                     if[SSAV( I ) < SSAV( I+1 ) ) RESULT( 18] = ULPINV;
                     IF[SSAV( I ) < ZERO ) RESULT( 18] = ULPINV;
                  } // 131
                  if ( MNMIN >= 1 ) {
                     if[SSAV( MNMIN ) < ZERO ) RESULT( 18] = ULPINV;
                  }
               }

               // Test CGEJSV
               // Note: CGEJSV only works for M >= N

               RESULT[19] = ZERO;
               RESULT[20] = ZERO;
               RESULT[21] = ZERO;
               RESULT[22] = ZERO;
               if ( M >= N ) {
                  IWTMP = 2*MNMIN*MNMIN + 2*MNMIN + max( M, N );
                  LSWORK = IWTMP + ( IWSPC-1 )*( LWORK-IWTMP ) / 3;
                  LSWORK = min( LSWORK, LWORK );
                  LSWORK = max( LSWORK, 1 );
                  if (IWSPC == 4) LSWORK = LWORK;
                  LRWORK = max( 7, N + 2*M);

                  clacpy('F', M, N, ASAV, LDA, VTSAV, LDA );
                 srnamc.SRNAMT = 'CGEJSV';
                  cgejsv('G', 'U', 'V', 'R', 'N', 'N', M, N, VTSAV, LDA, SSAV, USAV, LDU, A, LDVT, WORK, LWORK, RWORK, LRWORK, IWORK, IINFO );

                  // CGEJSV returns V not VH

                  for (J = 1; J <= N; J++) { // 133
                     for (I = 1; I <= N; I++) { // 132
                        VTSAV[J][I] = CONJG (A(I,J));
  132                END DO;
  133             END DO;

                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9995 )'GEJSV', IINFO, M, N, JTYPE, LSWORK, IOLDSD;
                     INFO = ( IINFO ).abs();
                     return;
                  }

                  // Do tests 19--22

                  cbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RWORK, RESULT( 19 ) );
                  if ( M != 0 && N != 0 ) {
                     cunt01('Columns', M, M, USAV, LDU, WORK, LWORK, RWORK, RESULT( 20 ) );
                     cunt01('Rows', N, N, VTSAV, LDVT, WORK, LWORK, RWORK, RESULT( 21 ) );
                  }
                  RESULT[22] = ZERO;
                  for (I = 1; I <= MNMIN - 1; I++) { // 134
                     if[SSAV( I ) < SSAV( I+1 ) ) RESULT( 22] = ULPINV;
                     IF[SSAV( I ) < ZERO ) RESULT( 22] = ULPINV;
                  } // 134
                  if ( MNMIN >= 1 ) {
                     if[SSAV( MNMIN ) < ZERO ) RESULT( 22] = ULPINV;
                  }
               }

               // Test CGESVDX

               // Factorize A

               clacpy('F', M, N, ASAV, LDA, A, LDA );
              srnamc.SRNAMT = 'CGESVDX';
               cgesvdx('V', 'V', 'A', M, N, A, LDA, VL, VU, IL, IU, NS, SSAV, USAV, LDU, VTSAV, LDVT, WORK, LWORK, RWORK, IWORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9995 )'GESVDX', IINFO, M, N, JTYPE, LSWORK, IOLDSD;
                  INFO = ( IINFO ).abs();
                  return;
               }

               // Do tests 1--4

               RESULT[23] = ZERO;
               RESULT[24] = ZERO;
               RESULT[25] = ZERO;
               cbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RWORK, RESULT( 23 ) );
               if ( M != 0 && N != 0 ) {
                  cunt01('Columns', MNMIN, M, USAV, LDU, WORK, LWORK, RWORK, RESULT( 24 ) );
                  cunt01('Rows', MNMIN, N, VTSAV, LDVT, WORK, LWORK, RWORK, RESULT( 25 ) );
               }
               RESULT[26] = ZERO;
               for (I = 1; I <= MNMIN - 1; I++) { // 140
                  if[SSAV( I ) < SSAV( I+1 ) ) RESULT( 26] = ULPINV;
                  IF[SSAV( I ) < ZERO ) RESULT( 26] = ULPINV;
               } // 140
               if ( MNMIN >= 1 ) {
                  if[SSAV( MNMIN ) < ZERO ) RESULT( 26] = ULPINV;
               }

               // Do partial SVDs, comparing to SSAV, USAV, and VTSAV

               RESULT[27] = ZERO;
               RESULT[28] = ZERO;
               RESULT[29] = ZERO;
               for (IJU = 0; IJU <= 1; IJU++) { // 170
                  for (IJVT = 0; IJVT <= 1; IJVT++) { // 160
                     if( ( IJU == 0 && IJVT == 0 ) || ( IJU == 1 && IJVT == 1 ) ) GO TO 160;
                     JOBU = CJOBV( IJU+1 );
                     JOBVT = CJOBV( IJVT+1 );
                     RANGE = CJOBR( 1 );
                     clacpy('F', M, N, ASAV, LDA, A, LDA );
                    srnamc.SRNAMT = 'CGESVDX';
                     cgesvdx(JOBU, JOBVT, 'A', M, N, A, LDA, VL, VU, IL, IU, NS, SSAV, U, LDU, VT, LDVT, WORK, LWORK, RWORK, IWORK, IINFO );

                     // Compare U

                     DIF = ZERO;
                     if ( M > 0 && N > 0 ) {
                        if ( IJU == 1 ) {
                           cunt03('C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, RWORK, DIF, IINFO );
                        }
                     }
                     RESULT[27] = max( RESULT( 27 ), DIF );

                     // Compare VT

                     DIF = ZERO;
                     if ( M > 0 && N > 0 ) {
                        if ( IJVT == 1 ) {
                           cunt03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, RWORK, DIF, IINFO );
                        }
                     }
                     RESULT[28] = max( RESULT( 28 ), DIF );

                     // Compare S

                     DIF = ZERO;
                     DIV = max( double( MNMIN )*ULP*S( 1 ), SLAMCH( 'Safe minimum' ) );
                     for (I = 1; I <= MNMIN - 1; I++) { // 150
                        if( SSAV( I ) < SSAV( I+1 ) ) DIF = ULPINV;
                        IF( SSAV( I ) < ZERO ) DIF = ULPINV;
                        DIF = max( DIF, ABS( SSAV( I )-S( I ) ) / DIV );
                     } // 150
                     RESULT[29] = max( RESULT( 29 ), DIF );
                  } // 160
               } // 170

               // Do tests 8--10

               for (I = 1; I <= 4; I++) { // 180
                  ISEED2[I] = ISEED( I );
               } // 180
               if ( MNMIN <= 1 ) {
                  IL = 1;
                  IU = max( 1, MNMIN );
               } else {
                  IL = 1 + INT( ( MNMIN-1 )*SLARND( 1, ISEED2 ) );
                  IU = 1 + INT( ( MNMIN-1 )*SLARND( 1, ISEED2 ) );
                  if ( IU < IL ) {
                     ITEMP = IU;
                     IU = IL;
                     IL = ITEMP;
                  }
               }
               clacpy('F', M, N, ASAV, LDA, A, LDA );
              srnamc.SRNAMT = 'CGESVDX';
               cgesvdx('V', 'V', 'I', M, N, A, LDA, VL, VU, IL, IU, NSI, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, IWORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9995 )'GESVDX', IINFO, M, N, JTYPE, LSWORK, IOLDSD;
                  INFO = ( IINFO ).abs();
                  return;
               }

               RESULT[30] = ZERO;
               RESULT[31] = ZERO;
               RESULT[32] = ZERO;
               cbdt05(M, N, ASAV, LDA, S, NSI, U, LDU, VT, LDVT, WORK, RESULT( 30 ) );
               if ( M != 0 && N != 0 ) {
                  cunt01('Columns', M, NSI, U, LDU, WORK, LWORK, RWORK, RESULT( 31 ) );
                  cunt01('Rows', NSI, N, VT, LDVT, WORK, LWORK, RWORK, RESULT( 32 ) );
               }

               // Do tests 11--13

               if ( MNMIN > 0 && NSI > 1 ) {
                  if ( IL != 1 ) {
                     VU = SSAV( IL ) + max( HALF*ABS( SSAV( IL )-SSAV( IL-1 ) ), ULP*ANORM, TWO*RTUNFL );
                  } else {
                     VU = SSAV( 1 ) + max( HALF*ABS( SSAV( NS )-SSAV( 1 ) ), ULP*ANORM, TWO*RTUNFL );
                  }
                  if ( IU != NS ) {
                     VL = SSAV( IU ) - max( ULP*ANORM, TWO*RTUNFL, HALF*ABS( SSAV( IU+1 )-SSAV( IU ) ) );
                  } else {
                     VL = SSAV( NS ) - max( ULP*ANORM, TWO*RTUNFL, HALF*ABS( SSAV( NS )-SSAV( 1 ) ) );
                  }
                  VL = max( VL,ZERO );
                  VU = max( VU,ZERO );
                  if (VL >= VU) VU = max( VU*2, VU+VL+HALF );
               } else {
                  VL = ZERO;
                  VU = ONE;
               }
               clacpy('F', M, N, ASAV, LDA, A, LDA );
              srnamc.SRNAMT = 'CGESVDX';
               cgesvdx('V', 'V', 'V', M, N, A, LDA, VL, VU, IL, IU, NSV, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, IWORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9995 )'GESVDX', IINFO, M, N, JTYPE, LSWORK, IOLDSD;
                  INFO = ( IINFO ).abs();
                  return;
               }

               RESULT[33] = ZERO;
               RESULT[34] = ZERO;
               RESULT[35] = ZERO;
               cbdt05(M, N, ASAV, LDA, S, NSV, U, LDU, VT, LDVT, WORK, RESULT( 33 ) );
               if ( M != 0 && N != 0 ) {
                  cunt01('Columns', M, NSV, U, LDU, WORK, LWORK, RWORK, RESULT( 34 ) );
                  cunt01('Rows', NSV, N, VT, LDVT, WORK, LWORK, RWORK, RESULT( 35 ) );
               }

               // End of Loop -- Check for RESULT(j) > THRESH

               NTEST = 0;
               NFAIL = 0;
               for (J = 1; J <= 39; J++) { // 190
                  if( RESULT( J ) >= ZERO ) NTEST = NTEST + 1;
                  IF( RESULT( J ) >= THRESH ) NFAIL = NFAIL + 1;
               } // 190

               if (NFAIL > 0) NTESTF = NTESTF + 1;
               if ( NTESTF == 1 ) {
                  WRITE( NOUNIT, FMT = 9999 );
                  WRITE( NOUNIT, FMT = 9998 )THRESH;
                  NTESTF = 2;
               }

               for (J = 1; J <= 39; J++) { // 200
                  if ( RESULT( J ) >= THRESH ) {
                     WRITE( NOUNIT, FMT = 9997 )M, N, JTYPE, IWSPC, IOLDSD, J, RESULT( J );
                  }
               } // 200

               NERRS = NERRS + NFAIL;
               NTESTT = NTESTT + NTEST;

            } // 290

         } // 300
      } // 310

      // Summary

      alasvm('CBD', NOUNIT, NERRS, NTESTT, 0 );

 9999 FORMAT( ' SVD -- Complex Singular Value Decomposition Driver \n Matrix types (see CDRVBD for details):', / / ' 1 = Zero matrix\n 2 = Identity matrix\n 3 = Evenly spaced singular values near 1\n 4 = Evenly spaced singular values near underflow\n 5 = Evenly spaced singular values near overflow', / / ' Tests performed: ( A is dense, U and V are unitary,\n${' ' * 19} S is an array, and Upartial, VTpartial, and\n${' ' * 19} Spartial are partially computed U, VT and S),', / );
 9998 FORMAT( ' Tests performed with Test Threshold = ', F8.2, / ' CGESVD: \n 1 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) \n 2 = | I - U**T U | / ( M ulp ) \n 3 = | I - VT VT**T | / ( N ulp ) \n 4 = 0 if S contains min(M,N) nonnegative values in decreasing order, else 1/ulp\n 5 = | U - Upartial | / ( M ulp )\n 6 = | VT - VTpartial | / ( N ulp )\n 7 = | S - Spartial | / ( min(M,N) ulp |S| )\n CGESDD: \n 8 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) \n 9 = | I - U**T U | / ( M ulp ) \n10 = | I - VT VT**T | / ( N ulp ) \n11 = 0 if S contains min(M,N) nonnegative values in decreasing order, else 1/ulp\n12 = | U - Upartial | / ( M ulp )\n13 = | VT - VTpartial | / ( N ulp )\n14 = | S - Spartial | / ( min(M,N) ulp |S| )\n CGESVJ: ', / / '15 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) \n16 = | I - U**T U | / ( M ulp ) \n17 = | I - VT VT**T | / ( N ulp ) \n18 = 0 if S contains min(M,N) nonnegative values in decreasing order, else 1/ulp\n CGESJV: ', / / '19 = | A - U diag(S) VT | / ( |A| max(M,N) ulp )',/ '20 = | I - U**T U | / ( M ulp ) ',/ '21 = | I - VT VT**T | / ( N ulp ) ',/ '22 = 0 if S contains min(M,N) nonnegative values in',' decreasing order, else 1/ulp',/ ' CGESVDX(V,V,A): ', /  '23 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ',/ '24 = | I - U**T U | / ( M ulp ) ',/ '25 = | I - VT VT**T | / ( N ulp ) ',/ '26 = 0 if S contains min(M,N) nonnegative values in',' decreasing order, else 1/ulp',/ '27 = | U - Upartial | / ( M ulp )',/ '28 = | VT - VTpartial | / ( N ulp )',/ '29 = | S - Spartial | / ( min(M,N) ulp |S| )',/ ' CGESVDX(V,V,I): ',/ '30 = | U**T A VT**T - diag(S) | / ( |A| max(M,N) ulp )',/ '31 = | I - U**T U | / ( M ulp ) ',/ '32 = | I - VT VT**T | / ( N ulp ) ',/ ' CGESVDX(V,V,V) ',/ '33 = | U**T A VT**T - diag(S) | / ( |A| max(M,N) ulp )',/ '34 = | I - U**T U | / ( M ulp ) ',/ '35 = | I - VT VT**T | / ( N ulp ) ',' CGESVDQ(H,N,N,A,A',/ '36 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ',/ '37 = | I - U**T U | / ( M ulp ) ',/ '38 = | I - VT VT**T | / ( N ulp ) ',/ '39 = 0 if S contains min(M,N) nonnegative values in',' decreasing order, else 1/ulp',/ / );
 9997 FORMAT( ' M=${.i5}, N=${.i5}, type ${.i1}, IWS=${.i1}, seed=${i4(4, ',')}', ' test(${.i2})=${.g11_4};
 9996 FORMAT( ' CDRVBD: ${} returned INFO=${.i6}.\n${' ' * 9}M=${.i6}, N=${.i6}, JTYPE=${.i6}, ISEED=(${i5(3, ',')}', I5, ')' );
 9995 FORMAT( ' CDRVBD: ${} returned INFO=${.i6}.\n${' ' * 9}M=${.i6}, N=${.i6}, JTYPE=${.i6}, LSWORK=', I6, / 9X, 'ISEED=(${i5(3, ',')}', I5, ')' );

      return;
      }
