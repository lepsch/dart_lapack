      void sdrvbd(NSIZES, MM, NN, NTYPES, DOTYPE, final Array<int> ISEED, THRESH, final Matrix<double> A, final int LDA, final Matrix<double> U, final int LDU, final Matrix<double> VT, final int LDVT, ASAV, USAV, VTSAV, S, SSAV, E, final Array<double> WORK, final int LWORK, IWORK, NOUT, Box<int> INFO ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      int                INFO, LDA, LDU, LDVT, LWORK, NOUT, NSIZES, NTYPES;
      double               THRESH;
      bool               DOTYPE( * );
      int                ISEED( 4 ), IWORK( * ), MM( * ), NN( * );
      double               A( LDA, * ), ASAV( LDA, * ), E( * ), S( * ), SSAV( * ), U( LDU, * ), USAV( LDU, * ), VT( LDVT, * ), VTSAV( LDVT, * ), WORK( * );
      // ..

      double              ZERO, ONE, TWO, HALF;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0, HALF = 0.5 ;
      int                MAXTYP;
      const              MAXTYP = 5 ;
      bool               BADMM, BADNN;
      String             JOBQ, JOBU, JOBVT, RANGE;
      String             PATH;
      int                I, IINFO, IJQ, IJU, IJVT, IL,IU, IWS, IWTMP, ITEMP, J, JSIZE, JTYPE, LSWORK, M, MINWRK, MMAX, MNMAX, MNMIN, MTYPES, N, NFAIL, NMAX, NS, NSI, NSV, NTEST;
      double               ANORM, DIF, DIV, OVFL, RTUNFL, ULP, ULPINV, UNFL, VL, VU;
      // ..
      // .. Local Scalars for DGESVDQ ..
      int                LIWORK, LRWORK, NUMRANK;
      // ..
      // .. Local Arrays for DGESVDQ ..
      double               RWORK( 2 );
      String             CJOB( 4 ), CJOBR( 3 ), CJOBV( 2 );
      int                IOLDSD( 4 ), ISEED2( 4 );
      double               RESULT( 39 );
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH, SLARND;
      // EXTERNAL SLAMCH, SLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALASVM, SBDT01, SGEJSV, SGESDD, SGESVD, SGESVDQ, SGESVDX, SGESVJ, SLACPY, SLASET, SLATMS, SORT01, SORT03, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, REAL, INT, MAX, MIN
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String            srnamc.SRNAMT;
      int                INFOT, NUNIT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NUNIT, OK, LERR
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Data statements ..
      const CJOB = [ 'N', 'O', 'S', 'A' ];
      const CJOBR = [ 'A', 'V', 'I' ];
      const CJOBV = [ 'N', 'V' ];

      // Check for errors

      INFO = 0;
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
         MINWRK = max( MINWRK, max( 3*min( MM( J ), NN( J ) )+max( MM( J ), NN( J ) ), 5*min( MM( J ), NN( J )-4 ) )+2*min( MM( J ), NN( J ) )**2 );
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
         xerbla('SDRVBD', -INFO );
         return;
      }

      // Initialize constants

      PATH[1: 1] = 'Single precision';
      PATH[2: 3] = 'BD';
      NFAIL = 0;
      NTEST = 0;
      UNFL = SLAMCH( 'Safe minimum' );
      OVFL = ONE / UNFL;
      ULP = SLAMCH( 'Precision' );
      RTUNFL = sqrt( UNFL );
      ULPINV = ONE / ULP;
      INFOT = 0;

      // Loop over sizes, types

      for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) { // 240
         M = MM( JSIZE );
         N = NN( JSIZE );
         MNMIN = min( M, N );

         if ( NSIZES != 1 ) {
            MTYPES = min( MAXTYP, NTYPES );
         } else {
            MTYPES = min( MAXTYP+1, NTYPES );
         }

         for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) { // 230
            if( !DOTYPE( JTYPE ) ) GO TO 230;

            for (J = 1; J <= 4; J++) { // 20
               IOLDSD[J] = ISEED( J );
            } // 20

            // Compute "A"

            if (MTYPES > MAXTYP) GO TO 30;

            if ( JTYPE == 1 ) {

               // Zero matrix

               slaset('Full', M, N, ZERO, ZERO, A, LDA );

            } else if ( JTYPE == 2 ) {

               // Identity matrix

               slaset('Full', M, N, ZERO, ONE, A, LDA );

            } else {

               // (Scaled) random matrix

               if (JTYPE == 3) ANORM = ONE;
               if( JTYPE == 4 ) ANORM = UNFL / ULP;
               IF( JTYPE == 5 ) ANORM = OVFL*ULP;
               slatms(M, N, 'U', ISEED, 'N', S, 4, double( MNMIN ), ANORM, M-1, N-1, 'N', A, LDA, WORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUT, FMT = 9996 )'Generator', IINFO, M, N, JTYPE, IOLDSD;
                  INFO = ( IINFO ).abs();
                  return;
               }
            }

            } // 30
            slacpy('F', M, N, A, LDA, ASAV, LDA );

            // Do for minimal and adequate (for blocking) workspace

            for (IWS = 1; IWS <= 4; IWS++) { // 220

               for (J = 1; J <= 32; J++) { // 40
                  RESULT[J] = -ONE;
               } // 40

               // Test SGESVD: Factorize A

               IWTMP = max( 3*min( M, N )+max( M, N ), 5*min( M, N ) );
               LSWORK = IWTMP + ( IWS-1 )*( LWORK-IWTMP ) / 3;
               LSWORK = min( LSWORK, LWORK );
               LSWORK = max( LSWORK, 1 );
               if (IWS == 4) LSWORK = LWORK;

               if (IWS > 1) slacpy( 'F', M, N, ASAV, LDA, A, LDA );
              srnamc.SRNAMT = 'SGESVD';
               sgesvd('A', 'A', M, N, A, LDA, SSAV, USAV, LDU, VTSAV, LDVT, WORK, LSWORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUT, FMT = 9995 )'GESVD', IINFO, M, N, JTYPE, LSWORK, IOLDSD;
                  INFO = ( IINFO ).abs();
                  return;
               }

               // Do tests 1--4

               sbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RESULT( 1 ) );
               if ( M != 0 && N != 0 ) {
                  sort01('Columns', M, M, USAV, LDU, WORK, LWORK, RESULT( 2 ) );
                  sort01('Rows', N, N, VTSAV, LDVT, WORK, LWORK, RESULT( 3 ) );
               }
               RESULT[4] = ZERO;
               for (I = 1; I <= MNMIN - 1; I++) { // 50
                  if[SSAV( I ) < SSAV( I+1 ) ) RESULT( 4] = ULPINV;
                  IF[SSAV( I ) < ZERO ) RESULT( 4] = ULPINV;
               } // 50
               if ( MNMIN >= 1 ) {
                  if[SSAV( MNMIN ) < ZERO ) RESULT( 4] = ULPINV;
               }

               // Do partial SVDs, comparing to SSAV, USAV, and VTSAV

               RESULT[5] = ZERO;
               RESULT[6] = ZERO;
               RESULT[7] = ZERO;
               for (IJU = 0; IJU <= 3; IJU++) { // 80
                  for (IJVT = 0; IJVT <= 3; IJVT++) { // 70
                     if( ( IJU == 3 && IJVT == 3 ) || ( IJU == 1 && IJVT == 1 ) )GO TO 70;
                     JOBU = CJOB( IJU+1 );
                     JOBVT = CJOB( IJVT+1 );
                     slacpy('F', M, N, ASAV, LDA, A, LDA );
                    srnamc.SRNAMT = 'SGESVD';
                     sgesvd(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LSWORK, IINFO );

                     // Compare U

                     DIF = ZERO;
                     if ( M > 0 && N > 0 ) {
                        if ( IJU == 1 ) {
                           sort03('C', M, MNMIN, M, MNMIN, USAV, LDU, A, LDA, WORK, LWORK, DIF, IINFO );
                        } else if ( IJU == 2 ) {
                           sort03('C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, DIF, IINFO );
                        } else if ( IJU == 3 ) {
                           sort03('C', M, M, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, DIF, IINFO );
                        }
                     }
                     RESULT[5] = max( RESULT( 5 ), DIF );

                     // Compare VT

                     DIF = ZERO;
                     if ( M > 0 && N > 0 ) {
                        if ( IJVT == 1 ) {
                           sort03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, A, LDA, WORK, LWORK, DIF, IINFO );
                        } else if ( IJVT == 2 ) {
                           sort03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, DIF, IINFO );
                        } else if ( IJVT == 3 ) {
                           sort03('R', N, N, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, DIF, IINFO );
                        }
                     }
                     RESULT[6] = max( RESULT( 6 ), DIF );

                     // Compare S

                     DIF = ZERO;
                     DIV = max( MNMIN*ULP*S( 1 ), UNFL );
                     for (I = 1; I <= MNMIN - 1; I++) { // 60
                        if( SSAV( I ) < SSAV( I+1 ) ) DIF = ULPINV;
                        IF( SSAV( I ) < ZERO ) DIF = ULPINV;
                        DIF = max( DIF, ABS( SSAV( I )-S( I ) ) / DIV );
                     } // 60
                     RESULT[7] = max( RESULT( 7 ), DIF );
                  } // 70
               } // 80

               // Test SGESDD: Factorize A

               IWTMP = 5*MNMIN*MNMIN + 9*MNMIN + max( M, N );
               LSWORK = IWTMP + ( IWS-1 )*( LWORK-IWTMP ) / 3;
               LSWORK = min( LSWORK, LWORK );
               LSWORK = max( LSWORK, 1 );
               if (IWS == 4) LSWORK = LWORK;

               slacpy('F', M, N, ASAV, LDA, A, LDA );
              srnamc.SRNAMT = 'SGESDD';
               sgesdd('A', M, N, A, LDA, SSAV, USAV, LDU, VTSAV, LDVT, WORK, LSWORK, IWORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUT, FMT = 9995 )'GESDD', IINFO, M, N, JTYPE, LSWORK, IOLDSD;
                  INFO = ( IINFO ).abs();
                  return;
               }

               // Do tests 8--11

               sbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RESULT( 8 ) );
               if ( M != 0 && N != 0 ) {
                  sort01('Columns', M, M, USAV, LDU, WORK, LWORK, RESULT( 9 ) );
                  sort01('Rows', N, N, VTSAV, LDVT, WORK, LWORK, RESULT( 10 ) );
               }
               RESULT[11] = ZERO;
               for (I = 1; I <= MNMIN - 1; I++) { // 90
                  if[SSAV( I ) < SSAV( I+1 ) ) RESULT( 11] = ULPINV;
                  IF[SSAV( I ) < ZERO ) RESULT( 11] = ULPINV;
               } // 90
               if ( MNMIN >= 1 ) {
                  if[SSAV( MNMIN ) < ZERO ) RESULT( 11] = ULPINV;
               }

               // Do partial SVDs, comparing to SSAV, USAV, and VTSAV

               RESULT[12] = ZERO;
               RESULT[13] = ZERO;
               RESULT[14] = ZERO;
               for (IJQ = 0; IJQ <= 2; IJQ++) { // 110
                  JOBQ = CJOB( IJQ+1 );
                  slacpy('F', M, N, ASAV, LDA, A, LDA );
                 srnamc.SRNAMT = 'SGESDD';
                  sgesdd(JOBQ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LSWORK, IWORK, IINFO );

                  // Compare U

                  DIF = ZERO;
                  if ( M > 0 && N > 0 ) {
                     if ( IJQ == 1 ) {
                        if ( M >= N ) {
                           sort03('C', M, MNMIN, M, MNMIN, USAV, LDU, A, LDA, WORK, LWORK, DIF, INFO );
                        } else {
                           sort03('C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, DIF, INFO );
                        }
                     } else if ( IJQ == 2 ) {
                        sort03('C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, DIF, INFO );
                     }
                  }
                  RESULT[12] = max( RESULT( 12 ), DIF );

                  // Compare VT

                  DIF = ZERO;
                  if ( M > 0 && N > 0 ) {
                     if ( IJQ == 1 ) {
                        if ( M >= N ) {
                           sort03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, DIF, INFO );
                        } else {
                           sort03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, A, LDA, WORK, LWORK, DIF, INFO );
                        }
                     } else if ( IJQ == 2 ) {
                        sort03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, DIF, INFO );
                     }
                  }
                  RESULT[13] = max( RESULT( 13 ), DIF );

                  // Compare S

                  DIF = ZERO;
                  DIV = max( MNMIN*ULP*S( 1 ), UNFL );
                  for (I = 1; I <= MNMIN - 1; I++) { // 100
                     if( SSAV( I ) < SSAV( I+1 ) ) DIF = ULPINV;
                     IF( SSAV( I ) < ZERO ) DIF = ULPINV;
                     DIF = max( DIF, ABS( SSAV( I )-S( I ) ) / DIV );
                  } // 100
                  RESULT[14] = max( RESULT( 14 ), DIF );
               } // 110

               // Test SGESVDQ
               // Note: SGESVDQ only works for M >= N

               RESULT[36] = ZERO;
               RESULT[37] = ZERO;
               RESULT[38] = ZERO;
               RESULT[39] = ZERO;

               if ( M >= N ) {
                  IWTMP = 5*MNMIN*MNMIN + 9*MNMIN + max( M, N );
                  LSWORK = IWTMP + ( IWS-1 )*( LWORK-IWTMP ) / 3;
                  LSWORK = min( LSWORK, LWORK );
                  LSWORK = max( LSWORK, 1 );
                  if (IWS == 4) LSWORK = LWORK;

                  slacpy('F', M, N, ASAV, LDA, A, LDA );
                 srnamc.SRNAMT = 'SGESVDQ';

                  LRWORK = 2;
                  LIWORK = max( N, 1 );
                  sgesvdq('H', 'N', 'N', 'A', 'A',  M, N, A, LDA, SSAV, USAV, LDU, VTSAV, LDVT, NUMRANK, IWORK, LIWORK, WORK, LWORK, RWORK, LRWORK, IINFO );

                  if ( IINFO != 0 ) {
                     WRITE( NOUT, FMT = 9995 )'SGESVDQ', IINFO, M, N, JTYPE, LSWORK, IOLDSD;
                     INFO = ( IINFO ).abs();
                     return;
                  }

                  // Do tests 36--39

                  sbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RESULT( 36 ) );
                  if ( M != 0 && N != 0 ) {
                     sort01('Columns', M, M, USAV, LDU, WORK, LWORK, RESULT( 37 ) );
                     sort01('Rows', N, N, VTSAV, LDVT, WORK, LWORK, RESULT( 38 ) );
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

               // Test SGESVJ
               // Note: SGESVJ only works for M >= N

               RESULT[15] = ZERO;
               RESULT[16] = ZERO;
               RESULT[17] = ZERO;
               RESULT[18] = ZERO;

               if ( M >= N ) {
                  IWTMP = 5*MNMIN*MNMIN + 9*MNMIN + max( M, N );
                  LSWORK = IWTMP + ( IWS-1 )*( LWORK-IWTMP ) / 3;
                  LSWORK = min( LSWORK, LWORK );
                  LSWORK = max( LSWORK, 1 );
                  if (IWS == 4) LSWORK = LWORK;

                  slacpy('F', M, N, ASAV, LDA, USAV, LDA );
                 srnamc.SRNAMT = 'SGESVJ';
                  sgesvj('G', 'U', 'V', M, N, USAV, LDA, SSAV, 0, A, LDVT, WORK, LWORK, INFO );

                  // SGESVJ returns V not VT

                  for (J = 1; J <= N; J++) {
                     for (I = 1; I <= N; I++) {
                        VTSAV[J][I] = A(I,J);
                     }
                  }

                  if ( IINFO != 0 ) {
                     WRITE( NOUT, FMT = 9995 )'GESVJ', IINFO, M, N, JTYPE, LSWORK, IOLDSD;
                     INFO = ( IINFO ).abs();
                     return;
                  }

                  // Do tests 15--18

                  sbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RESULT( 15 ) );
                  if ( M != 0 && N != 0 ) {
                     sort01('Columns', M, M, USAV, LDU, WORK, LWORK, RESULT( 16 ) );
                     sort01('Rows', N, N, VTSAV, LDVT, WORK, LWORK, RESULT( 17 ) );
                  }
                  RESULT[18] = ZERO;
                  for (I = 1; I <= MNMIN - 1; I++) { // 120
                     if[SSAV( I ) < SSAV( I+1 ) ) RESULT( 18] = ULPINV;
                     IF[SSAV( I ) < ZERO ) RESULT( 18] = ULPINV;
                  } // 120
                  if ( MNMIN >= 1 ) {
                     if[SSAV( MNMIN ) < ZERO ) RESULT( 18] = ULPINV;
                  }
               }

               // Test SGEJSV
               // Note: SGEJSV only works for M >= N

               RESULT[19] = ZERO;
               RESULT[20] = ZERO;
               RESULT[21] = ZERO;
               RESULT[22] = ZERO;
               if ( M >= N ) {
                  IWTMP = 5*MNMIN*MNMIN + 9*MNMIN + max( M, N );
                  LSWORK = IWTMP + ( IWS-1 )*( LWORK-IWTMP ) / 3;
                  LSWORK = min( LSWORK, LWORK );
                  LSWORK = max( LSWORK, 1 );
                  if (IWS == 4) LSWORK = LWORK;

                  slacpy('F', M, N, ASAV, LDA, VTSAV, LDA );
                 srnamc.SRNAMT = 'SGEJSV';
                  sgejsv('G', 'U', 'V', 'R', 'N', 'N', M, N, VTSAV, LDA, SSAV, USAV, LDU, A, LDVT, WORK, LWORK, IWORK, INFO );

                  // SGEJSV returns V not VT

                  for (J = 1; J <= N; J++) { // 140
                     for (I = 1; I <= N; I++) { // 130
                        VTSAV[J][I] = A(I,J);
  130                END DO;
  140             END DO;

                  if ( IINFO != 0 ) {
                     WRITE( NOUT, FMT = 9995 )'GEJSV', IINFO, M, N, JTYPE, LSWORK, IOLDSD;
                     INFO = ( IINFO ).abs();
                     return;
                  }

                  // Do tests 19--22

                  sbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RESULT( 19 ) );
                  if ( M != 0 && N != 0 ) {
                     sort01('Columns', M, M, USAV, LDU, WORK, LWORK, RESULT( 20 ) );
                     sort01('Rows', N, N, VTSAV, LDVT, WORK, LWORK, RESULT( 21 ) );
                  }
                  RESULT[22] = ZERO;
                  for (I = 1; I <= MNMIN - 1; I++) { // 150
                     if[SSAV( I ) < SSAV( I+1 ) ) RESULT( 22] = ULPINV;
                     IF[SSAV( I ) < ZERO ) RESULT( 22] = ULPINV;
                  } // 150
                  if ( MNMIN >= 1 ) {
                     if[SSAV( MNMIN ) < ZERO ) RESULT( 22] = ULPINV;
                  }
               }

               // Test SGESVDX

               slacpy('F', M, N, ASAV, LDA, A, LDA );
               sgesvdx('V', 'V', 'A', M, N, A, LDA, VL, VU, IL, IU, NS, SSAV, USAV, LDU, VTSAV, LDVT, WORK, LWORK, IWORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUT, FMT = 9995 )'GESVDX', IINFO, M, N, JTYPE, LSWORK, IOLDSD;
                  INFO = ( IINFO ).abs();
                  return;
               }

               // Do tests 23--29

               RESULT[23] = ZERO;
               RESULT[24] = ZERO;
               RESULT[25] = ZERO;
               sbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RESULT( 23 ) );
               if ( M != 0 && N != 0 ) {
                  sort01('Columns', M, M, USAV, LDU, WORK, LWORK, RESULT( 24 ) );
                  sort01('Rows', N, N, VTSAV, LDVT, WORK, LWORK, RESULT( 25 ) );
               }
               RESULT[26] = ZERO;
               for (I = 1; I <= MNMIN - 1; I++) { // 160
                  if[SSAV( I ) < SSAV( I+1 ) ) RESULT( 26] = ULPINV;
                  IF[SSAV( I ) < ZERO ) RESULT( 26] = ULPINV;
               } // 160
               if ( MNMIN >= 1 ) {
                  if[SSAV( MNMIN ) < ZERO ) RESULT( 26] = ULPINV;
               }

               // Do partial SVDs, comparing to SSAV, USAV, and VTSAV

               RESULT[27] = ZERO;
               RESULT[28] = ZERO;
               RESULT[29] = ZERO;
               for (IJU = 0; IJU <= 1; IJU++) { // 180
                  for (IJVT = 0; IJVT <= 1; IJVT++) { // 170
                     if( ( IJU == 0 && IJVT == 0 ) || ( IJU == 1 && IJVT == 1 ) )GO TO 170;
                     JOBU = CJOBV( IJU+1 );
                     JOBVT = CJOBV( IJVT+1 );
                     RANGE = CJOBR( 1 );
                     slacpy('F', M, N, ASAV, LDA, A, LDA );
                     sgesvdx(JOBU, JOBVT, RANGE, M, N, A, LDA, VL, VU, IL, IU, NS, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, IINFO );

                     // Compare U

                     DIF = ZERO;
                     if ( M > 0 && N > 0 ) {
                        if ( IJU == 1 ) {
                           sort03('C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, DIF, IINFO );
                        }
                     }
                     RESULT[27] = max( RESULT( 27 ), DIF );

                     // Compare VT

                     DIF = ZERO;
                     if ( M > 0 && N > 0 ) {
                        if ( IJVT == 1 ) {
                           sort03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, DIF, IINFO );
                        }
                     }
                     RESULT[28] = max( RESULT( 28 ), DIF );

                     // Compare S

                     DIF = ZERO;
                     DIV = max( MNMIN*ULP*S( 1 ), UNFL );
                     for (I = 1; I <= MNMIN - 1; I++) { // 190
                        if( SSAV( I ) < SSAV( I+1 ) ) DIF = ULPINV;
                        IF( SSAV( I ) < ZERO ) DIF = ULPINV;
                        DIF = max( DIF, ABS( SSAV( I )-S( I ) ) / DIV );
                     } // 190
                     RESULT[29] = max( RESULT( 29 ), DIF );
                  } // 170
               } // 180

               // Do tests 30--32: SGESVDX( 'V', 'V', 'I' )

               for (I = 1; I <= 4; I++) { // 200
                  ISEED2[I] = ISEED( I );
               } // 200
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
               slacpy('F', M, N, ASAV, LDA, A, LDA );
               sgesvdx('V', 'V', 'I', M, N, A, LDA, VL, VU, IL, IU, NSI, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUT, FMT = 9995 )'GESVDX', IINFO, M, N, JTYPE, LSWORK, IOLDSD;
                  INFO = ( IINFO ).abs();
                  return;
               }

               RESULT[30] = ZERO;
               RESULT[31] = ZERO;
               RESULT[32] = ZERO;
               sbdt05(M, N, ASAV, LDA, S, NSI, U, LDU, VT, LDVT, WORK, RESULT( 30 ) );
               sort01('Columns', M, NSI, U, LDU, WORK, LWORK, RESULT( 31 ) );
               sort01('Rows', NSI, N, VT, LDVT, WORK, LWORK, RESULT( 32 ) );

               // Do tests 33--35: SGESVDX( 'V', 'V', 'V' )

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
               slacpy('F', M, N, ASAV, LDA, A, LDA );
               sgesvdx('V', 'V', 'V', M, N, A, LDA, VL, VU, IL, IU, NSV, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUT, FMT = 9995 )'GESVDX', IINFO, M, N, JTYPE, LSWORK, IOLDSD;
                  INFO = ( IINFO ).abs();
                  return;
               }

               RESULT[33] = ZERO;
               RESULT[34] = ZERO;
               RESULT[35] = ZERO;
               sbdt05(M, N, ASAV, LDA, S, NSV, U, LDU, VT, LDVT, WORK, RESULT( 33 ) );
               sort01('Columns', M, NSV, U, LDU, WORK, LWORK, RESULT( 34 ) );
               sort01('Rows', NSV, N, VT, LDVT, WORK, LWORK, RESULT( 35 ) );

               // End of Loop -- Check for RESULT(j) > THRESH

               for (J = 1; J <= 39; J++) { // 210
                  if ( RESULT( J ) >= THRESH ) {
                     if ( NFAIL == 0 ) {
                        WRITE( NOUT, FMT = 9999 );
                        WRITE( NOUT, FMT = 9998 );
                     }
                     WRITE( NOUT, FMT = 9997 )M, N, JTYPE, IWS, IOLDSD, J, RESULT( J );
                     NFAIL = NFAIL + 1;
                  }
               } // 210
               NTEST = NTEST + 39;
            } // 220
         } // 230
      } // 240

      // Summary

      alasvm(PATH, NOUT, NFAIL, NTEST, 0 );

 9999 FORMAT( ' SVD -- Real Singular Value Decomposition Driver \n Matrix types (see SDRVBD for details):\n\n 1 = Zero matrix\n 2 = Identity matrix\n 3 = Evenly spaced singular values near 1\n 4 = Evenly spaced singular values near underflow\n 5 = Evenly spaced singular values near overflow\n\n Tests performed: ( A is dense, U and V are orthogonal,\n${' ' * 19} S is an array, and Upartial, VTpartial, and\n${' ' * 19} Spartial are partially computed U, VT and S),\n');
 9998 FORMAT( ' 1 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) \n 2 = | I - U**T U | / ( M ulp ) \n 3 = | I - VT VT**T | / ( N ulp ) \n 4 = 0 if S contains min(M,N) nonnegative values in decreasing order, else 1/ulp\n 5 = | U - Upartial | / ( M ulp )\n 6 = | VT - VTpartial | / ( N ulp )\n 7 = | S - Spartial | / ( min(M,N) ulp |S| )\n 8 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) \n 9 = | I - U**T U | / ( M ulp ) \n10 = | I - VT VT**T | / ( N ulp ) \n11 = 0 if S contains min(M,N) nonnegative values in decreasing order, else 1/ulp\n12 = | U - Upartial | / ( M ulp )\n13 = | VT - VTpartial | / ( N ulp )\n14 = | S - Spartial | / ( min(M,N) ulp |S| )\n15 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) \n16 = | I - U**T U | / ( M ulp ) \n17 = | I - VT VT**T | / ( N ulp ) \n18 = 0 if S contains min(M,N) nonnegative values in decreasing order, else 1/ulp\n19 = | U - Upartial | / ( M ulp )\n20 = | VT - VTpartial | / ( N ulp )\n21 = | S - Spartial | / ( min(M,N) ulp |S| )\n22 = 0 if S contains min(M,N) nonnegative values in',' decreasing order, else 1/ulp',' SGESVDX(V,V,A) \n23 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ),\n24 = | I - U**T U | / ( M ulp ) \n25 = | I - VT VT**T | / ( N ulp ) \n26 = 0 if S contains min(M,N) nonnegative values in',' decreasing order, else 1/ulp\n27 = | U - Upartial | / ( M ulp )\n28 = | VT - VTpartial | / ( N ulp )\n29 = | S - Spartial | / ( min(M,N) ulp |S| )\n30 = | U**T A VT**T - diag(S) | / ( |A| max(M,N) ulp ),',' SGESVDX(V,V,I) \n31 = | I - U**T U | / ( M ulp ) \n32 = | I - VT VT**T | / ( N ulp ) \n33 = | U**T A VT**T - diag(S) | / ( |A| max(M,N) ulp ),',' SGESVDX(V,V,V) \n34 = | I - U**T U | / ( M ulp ) \n35 = | I - VT VT**T | / ( N ulp ) ',' SGESVDQ(H,N,N,A,A\n36 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) \n37 = | I - U**T U | / ( M ulp ) \n38 = | I - VT VT**T | / ( N ulp ) \n39 = 0 if S contains min(M,N) nonnegative values in',' decreasing order, else 1/ulp',/ / );
 9997 FORMAT( ' M=${.i5}, N=${.i5}, type ${.i1}, IWS=${.i1}, seed=${i4(4, ',')}', ' test(${.i2})=${.g11_4};
 9996 FORMAT( ' SDRVBD: ${} returned INFO=${.i6}.\n${' ' * 9}M=${.i6}, N=${.i6}, JTYPE=${.i6}, ISEED=(${.i5(4, ',')})' );
 9995 FORMAT( ' SDRVBD: ${} returned INFO=${.i6}.\n${' ' * 9}M=${.i6}, N=${.i6}, JTYPE=${.i6}, LSWORK=', I6, / 9X, 'ISEED=(${.i5(4, ',')})' );
      }
