      SUBROUTINE DDRVBD( NSIZES, MM, NN, NTYPES, DOTYPE, ISEED, THRESH, A, LDA, U, LDU, VT, LDVT, ASAV, USAV, VTSAV, S, SSAV, E, WORK, LWORK, IWORK, NOUT, INFO )

      IMPLICIT NONE

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDU, LDVT, LWORK, NOUT, NSIZES, NTYPES;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                ISEED( 4 ), IWORK( * ), MM( * ), NN( * );
      double             A( LDA, * ), ASAV( LDA, * ), E( * ), S( * ), SSAV( * ), U( LDU, * ), USAV( LDU, * ), VT( LDVT, * ), VTSAV( LDVT, * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double            ZERO, ONE, TWO, HALF;
      const              ZERO = 0.0, ONE = 1.0, TWO = 2.0, HALF = 0.5 ;
      int                MAXTYP;
      const              MAXTYP = 5 ;
      // ..
      // .. Local Scalars ..
      bool               BADMM, BADNN;
      String             JOBQ, JOBU, JOBVT, RANGE;
      String             PATH;
      int                I, IINFO, IJQ, IJU, IJVT, IL,IU, IWS, IWTMP, ITEMP, J, JSIZE, JTYPE, LSWORK, M, MINWRK, MMAX, MNMAX, MNMIN, MTYPES, N, NFAIL, NMAX, NS, NSI, NSV, NTEST;
      double             ANORM, DIF, DIV, OVFL, RTUNFL, ULP, ULPINV, UNFL, VL, VU;
      // ..
      // .. Local Scalars for DGESVDQ ..
      int                LIWORK, LRWORK, NUMRANK;
      // ..
      // .. Local Arrays for DGESVDQ ..
      double             RWORK( 2 );
      // ..
      // .. Local Arrays ..
      String             CJOB( 4 ), CJOBR( 3 ), CJOBV( 2 );
      int                IOLDSD( 4 ), ISEED2( 4 );
      double             RESULT( 39 );
      // ..
      // .. External Functions ..
      double             DLAMCH, DLARND;
      // EXTERNAL DLAMCH, DLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALASVM, DBDT01, DGEJSV, DGESDD, DGESVD, DGESVDQ, DGESVDX, DGESVJ, DLACPY, DLASET, DLATMS, DORT01, DORT03, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, INT, MAX, MIN
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NUNIT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NUNIT, OK, LERR
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Data statements ..
      DATA               CJOB / 'N', 'O', 'S', 'A' /
      DATA               CJOBR / 'A', 'V', 'I' /
      DATA               CJOBV / 'N', 'V' /
      // ..
      // .. Executable Statements ..

      // Check for errors

      INFO = 0
      BADMM = false;
      BADNN = false;
      MMAX = 1
      NMAX = 1
      MNMAX = 1
      MINWRK = 1
      for (J = 1; J <= NSIZES; J++) { // 10
         MMAX = MAX( MMAX, MM( J ) )
         IF( MM( J ) < 0 ) BADMM = true;
         NMAX = MAX( NMAX, NN( J ) )
         IF( NN( J ) < 0 ) BADNN = true;
         MNMAX = MAX( MNMAX, MIN( MM( J ), NN( J ) ) )
         MINWRK = MAX( MINWRK, MAX( 3*MIN( MM( J ), NN( J ) )+MAX( MM( J ), NN( J ) ), 5*MIN( MM( J ), NN( J )-4 ) )+2*MIN( MM( J ), NN( J ) )**2 )
      } // 10

      // Check for errors

      if ( NSIZES < 0 ) {
         INFO = -1
      } else if ( BADMM ) {
         INFO = -2
      } else if ( BADNN ) {
         INFO = -3
      } else if ( NTYPES < 0 ) {
         INFO = -4
      } else if ( LDA < MAX( 1, MMAX ) ) {
         INFO = -10
      } else if ( LDU < MAX( 1, MMAX ) ) {
         INFO = -12
      } else if ( LDVT < MAX( 1, NMAX ) ) {
         INFO = -14
      } else if ( MINWRK > LWORK ) {
         INFO = -21
      }

      if ( INFO != 0 ) {
         xerbla('DDRVBD', -INFO );
         RETURN
      }

      // Initialize constants

      PATH( 1: 1 ) = 'double          ';
      PATH( 2: 3 ) = 'BD'
      NFAIL = 0
      NTEST = 0
      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      ULP = DLAMCH( 'Precision' )
      RTUNFL = SQRT( UNFL )
      ULPINV = ONE / ULP
      INFOT = 0

      // Loop over sizes, types

      for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) { // 240
         M = MM( JSIZE )
         N = NN( JSIZE )
         MNMIN = MIN( M, N )

         if ( NSIZES != 1 ) {
            MTYPES = MIN( MAXTYP, NTYPES )
         } else {
            MTYPES = MIN( MAXTYP+1, NTYPES )
         }

         for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) { // 230
            IF( .NOT.DOTYPE( JTYPE ) ) GO TO 230

            for (J = 1; J <= 4; J++) { // 20
               IOLDSD( J ) = ISEED( J )
            } // 20

            // Compute "A"

            if (MTYPES > MAXTYP) GO TO 30;

            if ( JTYPE == 1 ) {

               // Zero matrix

               dlaset('Full', M, N, ZERO, ZERO, A, LDA );

            } else if ( JTYPE == 2 ) {

               // Identity matrix

               dlaset('Full', M, N, ZERO, ONE, A, LDA );

            } else {

               // (Scaled) random matrix

               if (JTYPE == 3) ANORM = ONE                IF( JTYPE == 4 ) ANORM = UNFL / ULP                IF( JTYPE == 5 ) ANORM = OVFL*ULP;
               dlatms(M, N, 'U', ISEED, 'N', S, 4, DBLE( MNMIN ), ANORM, M-1, N-1, 'N', A, LDA, WORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUT, FMT = 9996 )'Generator', IINFO, M, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               }
            }

            } // 30
            dlacpy('F', M, N, A, LDA, ASAV, LDA );

            // Do for minimal and adequate (for blocking) workspace

            for (IWS = 1; IWS <= 4; IWS++) { // 220

               for (J = 1; J <= 32; J++) { // 40
                  RESULT( J ) = -ONE
               } // 40

               // Test DGESVD: Factorize A

               IWTMP = MAX( 3*MIN( M, N )+MAX( M, N ), 5*MIN( M, N ) )
               LSWORK = IWTMP + ( IWS-1 )*( LWORK-IWTMP ) / 3
               LSWORK = MIN( LSWORK, LWORK )
               LSWORK = MAX( LSWORK, 1 )
               if (IWS == 4) LSWORK = LWORK;

               if (IWS > 1) CALL DLACPY( 'F', M, N, ASAV, LDA, A, LDA );
               SRNAMT = 'DGESVD'
               dgesvd('A', 'A', M, N, A, LDA, SSAV, USAV, LDU, VTSAV, LDVT, WORK, LSWORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUT, FMT = 9995 )'GESVD', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               }

               // Do tests 1--4

               dbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RESULT( 1 ) );
               if ( M != 0 && N != 0 ) {
                  dort01('Columns', M, M, USAV, LDU, WORK, LWORK, RESULT( 2 ) );
                  dort01('Rows', N, N, VTSAV, LDVT, WORK, LWORK, RESULT( 3 ) );
               }
               RESULT( 4 ) = ZERO
               for (I = 1; I <= MNMIN - 1; I++) { // 50
                  IF( SSAV( I ) < SSAV( I+1 ) ) RESULT( 4 ) = ULPINV                   IF( SSAV( I ) < ZERO ) RESULT( 4 ) = ULPINV
               } // 50
               if ( MNMIN >= 1 ) {
                  IF( SSAV( MNMIN ) < ZERO ) RESULT( 4 ) = ULPINV
               }

               // Do partial SVDs, comparing to SSAV, USAV, and VTSAV

               RESULT( 5 ) = ZERO
               RESULT( 6 ) = ZERO
               RESULT( 7 ) = ZERO
               for (IJU = 0; IJU <= 3; IJU++) { // 80
                  for (IJVT = 0; IJVT <= 3; IJVT++) { // 70
                     IF( ( IJU == 3 && IJVT == 3 ) || ( IJU == 1 && IJVT == 1 ) )GO TO 70
                     JOBU = CJOB( IJU+1 )
                     JOBVT = CJOB( IJVT+1 )
                     dlacpy('F', M, N, ASAV, LDA, A, LDA );
                     SRNAMT = 'DGESVD'
                     dgesvd(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LSWORK, IINFO );

                     // Compare U

                     DIF = ZERO
                     if ( M > 0 && N > 0 ) {
                        if ( IJU == 1 ) {
                           dort03('C', M, MNMIN, M, MNMIN, USAV, LDU, A, LDA, WORK, LWORK, DIF, IINFO );
                        } else if ( IJU == 2 ) {
                           dort03('C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, DIF, IINFO );
                        } else if ( IJU == 3 ) {
                           dort03('C', M, M, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, DIF, IINFO );
                        }
                     }
                     RESULT( 5 ) = MAX( RESULT( 5 ), DIF )

                     // Compare VT

                     DIF = ZERO
                     if ( M > 0 && N > 0 ) {
                        if ( IJVT == 1 ) {
                           dort03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, A, LDA, WORK, LWORK, DIF, IINFO );
                        } else if ( IJVT == 2 ) {
                           dort03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, DIF, IINFO );
                        } else if ( IJVT == 3 ) {
                           dort03('R', N, N, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, DIF, IINFO );
                        }
                     }
                     RESULT( 6 ) = MAX( RESULT( 6 ), DIF )

                     // Compare S

                     DIF = ZERO
                     DIV = MAX( MNMIN*ULP*S( 1 ), UNFL )
                     for (I = 1; I <= MNMIN - 1; I++) { // 60
                        IF( SSAV( I ) < SSAV( I+1 ) ) DIF = ULPINV                         IF( SSAV( I ) < ZERO ) DIF = ULPINV
                        DIF = MAX( DIF, ABS( SSAV( I )-S( I ) ) / DIV )
                     } // 60
                     RESULT( 7 ) = MAX( RESULT( 7 ), DIF )
                  } // 70
               } // 80

               // Test DGESDD: Factorize A

               IWTMP = 5*MNMIN*MNMIN + 9*MNMIN + MAX( M, N )
               LSWORK = IWTMP + ( IWS-1 )*( LWORK-IWTMP ) / 3
               LSWORK = MIN( LSWORK, LWORK )
               LSWORK = MAX( LSWORK, 1 )
               if (IWS == 4) LSWORK = LWORK;

               dlacpy('F', M, N, ASAV, LDA, A, LDA );
               SRNAMT = 'DGESDD'
               dgesdd('A', M, N, A, LDA, SSAV, USAV, LDU, VTSAV, LDVT, WORK, LSWORK, IWORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUT, FMT = 9995 )'GESDD', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               }

               // Do tests 8--11

               dbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RESULT( 8 ) );
               if ( M != 0 && N != 0 ) {
                  dort01('Columns', M, M, USAV, LDU, WORK, LWORK, RESULT( 9 ) );
                  dort01('Rows', N, N, VTSAV, LDVT, WORK, LWORK, RESULT( 10 ) );
               }
               RESULT( 11 ) = ZERO
               for (I = 1; I <= MNMIN - 1; I++) { // 90
                  IF( SSAV( I ) < SSAV( I+1 ) ) RESULT( 11 ) = ULPINV                   IF( SSAV( I ) < ZERO ) RESULT( 11 ) = ULPINV
               } // 90
               if ( MNMIN >= 1 ) {
                  IF( SSAV( MNMIN ) < ZERO ) RESULT( 11 ) = ULPINV
               }

               // Do partial SVDs, comparing to SSAV, USAV, and VTSAV

               RESULT( 12 ) = ZERO
               RESULT( 13 ) = ZERO
               RESULT( 14 ) = ZERO
               for (IJQ = 0; IJQ <= 2; IJQ++) { // 110
                  JOBQ = CJOB( IJQ+1 )
                  dlacpy('F', M, N, ASAV, LDA, A, LDA );
                  SRNAMT = 'DGESDD'
                  dgesdd(JOBQ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LSWORK, IWORK, IINFO );

                  // Compare U

                  DIF = ZERO
                  if ( M > 0 && N > 0 ) {
                     if ( IJQ == 1 ) {
                        if ( M >= N ) {
                           dort03('C', M, MNMIN, M, MNMIN, USAV, LDU, A, LDA, WORK, LWORK, DIF, INFO );
                        } else {
                           dort03('C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, DIF, INFO );
                        }
                     } else if ( IJQ == 2 ) {
                        dort03('C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, DIF, INFO );
                     }
                  }
                  RESULT( 12 ) = MAX( RESULT( 12 ), DIF )

                  // Compare VT

                  DIF = ZERO
                  if ( M > 0 && N > 0 ) {
                     if ( IJQ == 1 ) {
                        if ( M >= N ) {
                           dort03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, DIF, INFO );
                        } else {
                           dort03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, A, LDA, WORK, LWORK, DIF, INFO );
                        }
                     } else if ( IJQ == 2 ) {
                        dort03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, DIF, INFO );
                     }
                  }
                  RESULT( 13 ) = MAX( RESULT( 13 ), DIF )

                  // Compare S

                  DIF = ZERO
                  DIV = MAX( MNMIN*ULP*S( 1 ), UNFL )
                  for (I = 1; I <= MNMIN - 1; I++) { // 100
                     IF( SSAV( I ) < SSAV( I+1 ) ) DIF = ULPINV                      IF( SSAV( I ) < ZERO ) DIF = ULPINV
                     DIF = MAX( DIF, ABS( SSAV( I )-S( I ) ) / DIV )
                  } // 100
                  RESULT( 14 ) = MAX( RESULT( 14 ), DIF )
               } // 110

               // Test DGESVDQ
               // Note: DGESVDQ only works for M >= N

               RESULT( 36 ) = ZERO
               RESULT( 37 ) = ZERO
               RESULT( 38 ) = ZERO
               RESULT( 39 ) = ZERO

               if ( M >= N ) {
                  IWTMP = 5*MNMIN*MNMIN + 9*MNMIN + MAX( M, N )
                  LSWORK = IWTMP + ( IWS-1 )*( LWORK-IWTMP ) / 3
                  LSWORK = MIN( LSWORK, LWORK )
                  LSWORK = MAX( LSWORK, 1 )
                  if (IWS == 4) LSWORK = LWORK;

                  dlacpy('F', M, N, ASAV, LDA, A, LDA );
                  SRNAMT = 'DGESVDQ'

                  LRWORK = 2
                  LIWORK = MAX( N, 1 )
                  dgesvdq('H', 'N', 'N', 'A', 'A',  M, N, A, LDA, SSAV, USAV, LDU, VTSAV, LDVT, NUMRANK, IWORK, LIWORK, WORK, LWORK, RWORK, LRWORK, IINFO );

                  if ( IINFO != 0 ) {
                     WRITE( NOUT, FMT = 9995 )'DGESVDQ', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                     INFO = ABS( IINFO )
                     RETURN
                  }

                  // Do tests 36--39

                  dbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RESULT( 36 ) );
                  if ( M != 0 && N != 0 ) {
                     dort01('Columns', M, M, USAV, LDU, WORK, LWORK, RESULT( 37 ) );
                     dort01('Rows', N, N, VTSAV, LDVT, WORK, LWORK, RESULT( 38 ) );
                  }
                  RESULT( 39 ) = ZERO
                  for (I = 1; I <= MNMIN - 1; I++) { // 199
                     IF( SSAV( I ) < SSAV( I+1 ) ) RESULT( 39 ) = ULPINV                      IF( SSAV( I ) < ZERO ) RESULT( 39 ) = ULPINV
                  } // 199
                  if ( MNMIN >= 1 ) {
                     IF( SSAV( MNMIN ) < ZERO ) RESULT( 39 ) = ULPINV
                  }
               }

               // Test DGESVJ
               // Note: DGESVJ only works for M >= N

               RESULT( 15 ) = ZERO
               RESULT( 16 ) = ZERO
               RESULT( 17 ) = ZERO
               RESULT( 18 ) = ZERO

               if ( M >= N ) {
                  IWTMP = 5*MNMIN*MNMIN + 9*MNMIN + MAX( M, N )
                  LSWORK = IWTMP + ( IWS-1 )*( LWORK-IWTMP ) / 3
                  LSWORK = MIN( LSWORK, LWORK )
                  LSWORK = MAX( LSWORK, 1 )
                  if (IWS == 4) LSWORK = LWORK;

                  dlacpy('F', M, N, ASAV, LDA, USAV, LDA );
                  SRNAMT = 'DGESVJ'
                  dgesvj('G', 'U', 'V', M, N, USAV, LDA, SSAV, 0, A, LDVT, WORK, LWORK, INFO );

                  // DGESVJ returns V not VT

                  for (J = 1; J <= N; J++) {
                     for (I = 1; I <= N; I++) {
                        VTSAV(J,I) = A(I,J)
                     }
                  }

                  if ( IINFO != 0 ) {
                     WRITE( NOUT, FMT = 9995 )'GESVJ', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                     INFO = ABS( IINFO )
                     RETURN
                  }

                  // Do tests 15--18

                  dbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RESULT( 15 ) );
                  if ( M != 0 && N != 0 ) {
                     dort01('Columns', M, M, USAV, LDU, WORK, LWORK, RESULT( 16 ) );
                     dort01('Rows', N, N, VTSAV, LDVT, WORK, LWORK, RESULT( 17 ) );
                  }
                  RESULT( 18 ) = ZERO
                  for (I = 1; I <= MNMIN - 1; I++) { // 120
                     IF( SSAV( I ) < SSAV( I+1 ) ) RESULT( 18 ) = ULPINV                      IF( SSAV( I ) < ZERO ) RESULT( 18 ) = ULPINV
                  } // 120
                  if ( MNMIN >= 1 ) {
                     IF( SSAV( MNMIN ) < ZERO ) RESULT( 18 ) = ULPINV
                  }
               }

               // Test DGEJSV
               // Note: DGEJSV only works for M >= N

               RESULT( 19 ) = ZERO
               RESULT( 20 ) = ZERO
               RESULT( 21 ) = ZERO
               RESULT( 22 ) = ZERO
               if ( M >= N ) {
                  IWTMP = 5*MNMIN*MNMIN + 9*MNMIN + MAX( M, N )
                  LSWORK = IWTMP + ( IWS-1 )*( LWORK-IWTMP ) / 3
                  LSWORK = MIN( LSWORK, LWORK )
                  LSWORK = MAX( LSWORK, 1 )
                  if (IWS == 4) LSWORK = LWORK;

                  dlacpy('F', M, N, ASAV, LDA, VTSAV, LDA );
                  SRNAMT = 'DGEJSV'
                  dgejsv('G', 'U', 'V', 'R', 'N', 'N', M, N, VTSAV, LDA, SSAV, USAV, LDU, A, LDVT, WORK, LWORK, IWORK, INFO );

                  // DGEJSV returns V not VT

                  for (J = 1; J <= N; J++) { // 140
                     for (I = 1; I <= N; I++) { // 130
                        VTSAV(J,I) = A(I,J)
  130                END DO
  140             END DO

                  if ( IINFO != 0 ) {
                     WRITE( NOUT, FMT = 9995 )'GEJSV', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                     INFO = ABS( IINFO )
                     RETURN
                  }

                  // Do tests 19--22

                  dbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RESULT( 19 ) );
                  if ( M != 0 && N != 0 ) {
                     dort01('Columns', M, M, USAV, LDU, WORK, LWORK, RESULT( 20 ) );
                     dort01('Rows', N, N, VTSAV, LDVT, WORK, LWORK, RESULT( 21 ) );
                  }
                  RESULT( 22 ) = ZERO
                  for (I = 1; I <= MNMIN - 1; I++) { // 150
                     IF( SSAV( I ) < SSAV( I+1 ) ) RESULT( 22 ) = ULPINV                      IF( SSAV( I ) < ZERO ) RESULT( 22 ) = ULPINV
                  } // 150
                  if ( MNMIN >= 1 ) {
                     IF( SSAV( MNMIN ) < ZERO ) RESULT( 22 ) = ULPINV
                  }
               }

               // Test DGESVDX

               dlacpy('F', M, N, ASAV, LDA, A, LDA );
               dgesvdx('V', 'V', 'A', M, N, A, LDA, VL, VU, IL, IU, NS, SSAV, USAV, LDU, VTSAV, LDVT, WORK, LWORK, IWORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUT, FMT = 9995 )'GESVDX', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               }

               // Do tests 23--29

               RESULT( 23 ) = ZERO
               RESULT( 24 ) = ZERO
               RESULT( 25 ) = ZERO
               dbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RESULT( 23 ) );
               if ( M != 0 && N != 0 ) {
                  dort01('Columns', M, M, USAV, LDU, WORK, LWORK, RESULT( 24 ) );
                  dort01('Rows', N, N, VTSAV, LDVT, WORK, LWORK, RESULT( 25 ) );
               }
               RESULT( 26 ) = ZERO
               for (I = 1; I <= MNMIN - 1; I++) { // 160
                  IF( SSAV( I ) < SSAV( I+1 ) ) RESULT( 26 ) = ULPINV                   IF( SSAV( I ) < ZERO ) RESULT( 26 ) = ULPINV
               } // 160
               if ( MNMIN >= 1 ) {
                  IF( SSAV( MNMIN ) < ZERO ) RESULT( 26 ) = ULPINV
               }

               // Do partial SVDs, comparing to SSAV, USAV, and VTSAV

               RESULT( 27 ) = ZERO
               RESULT( 28 ) = ZERO
               RESULT( 29 ) = ZERO
               for (IJU = 0; IJU <= 1; IJU++) { // 180
                  for (IJVT = 0; IJVT <= 1; IJVT++) { // 170
                     IF( ( IJU == 0 && IJVT == 0 ) || ( IJU == 1 && IJVT == 1 ) )GO TO 170
                     JOBU = CJOBV( IJU+1 )
                     JOBVT = CJOBV( IJVT+1 )
                     RANGE = CJOBR( 1 )
                     dlacpy('F', M, N, ASAV, LDA, A, LDA );
                     dgesvdx(JOBU, JOBVT, RANGE, M, N, A, LDA, VL, VU, IL, IU, NS, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, IINFO );

                     // Compare U

                     DIF = ZERO
                     if ( M > 0 && N > 0 ) {
                        if ( IJU == 1 ) {
                           dort03('C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, DIF, IINFO );
                        }
                     }
                     RESULT( 27 ) = MAX( RESULT( 27 ), DIF )

                     // Compare VT

                     DIF = ZERO
                     if ( M > 0 && N > 0 ) {
                        if ( IJVT == 1 ) {
                           dort03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, DIF, IINFO );
                        }
                     }
                     RESULT( 28 ) = MAX( RESULT( 28 ), DIF )

                     // Compare S

                     DIF = ZERO
                     DIV = MAX( MNMIN*ULP*S( 1 ), UNFL )
                     for (I = 1; I <= MNMIN - 1; I++) { // 190
                        IF( SSAV( I ) < SSAV( I+1 ) ) DIF = ULPINV                         IF( SSAV( I ) < ZERO ) DIF = ULPINV
                        DIF = MAX( DIF, ABS( SSAV( I )-S( I ) ) / DIV )
                     } // 190
                     RESULT( 29 ) = MAX( RESULT( 29 ), DIF )
                  } // 170
               } // 180

               // Do tests 30--32: DGESVDX( 'V', 'V', 'I' )

               for (I = 1; I <= 4; I++) { // 200
                  ISEED2( I ) = ISEED( I )
               } // 200
               if ( MNMIN <= 1 ) {
                  IL = 1
                  IU = MAX( 1, MNMIN )
               } else {
                  IL = 1 + INT( ( MNMIN-1 )*DLARND( 1, ISEED2 ) )
                  IU = 1 + INT( ( MNMIN-1 )*DLARND( 1, ISEED2 ) )
                  if ( IU < IL ) {
                     ITEMP = IU
                     IU = IL
                     IL = ITEMP
                  }
               }
               dlacpy('F', M, N, ASAV, LDA, A, LDA );
               dgesvdx('V', 'V', 'I', M, N, A, LDA, VL, VU, IL, IU, NSI, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUT, FMT = 9995 )'GESVDX', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               }

               RESULT( 30 ) = ZERO
               RESULT( 31 ) = ZERO
               RESULT( 32 ) = ZERO
               dbdt05(M, N, ASAV, LDA, S, NSI, U, LDU, VT, LDVT, WORK, RESULT( 30 ) );
               dort01('Columns', M, NSI, U, LDU, WORK, LWORK, RESULT( 31 ) );
               dort01('Rows', NSI, N, VT, LDVT, WORK, LWORK, RESULT( 32 ) );

               // Do tests 33--35: DGESVDX( 'V', 'V', 'V' )

               if ( MNMIN > 0 && NSI > 1 ) {
                  if ( IL != 1 ) {
                     VU = SSAV( IL ) + MAX( HALF*ABS( SSAV( IL )-SSAV( IL-1 ) ), ULP*ANORM, TWO*RTUNFL )
                  } else {
                     VU = SSAV( 1 ) + MAX( HALF*ABS( SSAV( NS )-SSAV( 1 ) ), ULP*ANORM, TWO*RTUNFL )
                  }
                  if ( IU != NS ) {
                     VL = SSAV( IU ) - MAX( ULP*ANORM, TWO*RTUNFL, HALF*ABS( SSAV( IU+1 )-SSAV( IU ) ) )
                  } else {
                     VL = SSAV( NS ) - MAX( ULP*ANORM, TWO*RTUNFL, HALF*ABS( SSAV( NS )-SSAV( 1 ) ) )
                  }
                  VL = MAX( VL,ZERO )
                  VU = MAX( VU,ZERO )
                  if (VL >= VU) VU = MAX( VU*2, VU+VL+HALF );
               } else {
                  VL = ZERO
                  VU = ONE
               }
               dlacpy('F', M, N, ASAV, LDA, A, LDA );
               dgesvdx('V', 'V', 'V', M, N, A, LDA, VL, VU, IL, IU, NSV, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUT, FMT = 9995 )'GESVDX', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               }

               RESULT( 33 ) = ZERO
               RESULT( 34 ) = ZERO
               RESULT( 35 ) = ZERO
               dbdt05(M, N, ASAV, LDA, S, NSV, U, LDU, VT, LDVT, WORK, RESULT( 33 ) );
               dort01('Columns', M, NSV, U, LDU, WORK, LWORK, RESULT( 34 ) );
               dort01('Rows', NSV, N, VT, LDVT, WORK, LWORK, RESULT( 35 ) );

               // End of Loop -- Check for RESULT(j) > THRESH

               for (J = 1; J <= 39; J++) { // 210
                  if ( RESULT( J ) >= THRESH ) {
                     if ( NFAIL == 0 ) {
                        WRITE( NOUT, FMT = 9999 )
                        WRITE( NOUT, FMT = 9998 )
                     }
                     WRITE( NOUT, FMT = 9997 )M, N, JTYPE, IWS, IOLDSD, J, RESULT( J )
                     NFAIL = NFAIL + 1
                  }
               } // 210
               NTEST = NTEST + 39
            } // 220
         } // 230
      } // 240

      // Summary

      alasvm(PATH, NOUT, NFAIL, NTEST, 0 );

 9999 FORMAT( ' SVD -- Real Singular Value Decomposition Driver ', / ' Matrix types (see DDRVBD for details):', / / ' 1 = Zero matrix', / ' 2 = Identity matrix', / ' 3 = Evenly spaced singular values near 1', / ' 4 = Evenly spaced singular values near underflow', / ' 5 = Evenly spaced singular values near overflow', / / ' Tests performed: ( A is dense, U and V are orthogonal,', / 19X, ' S is an array, and Upartial, VTpartial, and', / 19X, ' Spartial are partially computed U, VT and S),', / )
 9998 FORMAT( ' 1 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ', / ' 2 = | I - U**T U | / ( M ulp ) ', / ' 3 = | I - VT VT**T | / ( N ulp ) ', / ' 4 = 0 if S contains min(M,N) nonnegative values in', ' decreasing order, else 1/ulp', / ' 5 = | U - Upartial | / ( M ulp )', / ' 6 = | VT - VTpartial | / ( N ulp )', / ' 7 = | S - Spartial | / ( min(M,N) ulp |S| )', / ' 8 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ', / ' 9 = | I - U**T U | / ( M ulp ) ', / '10 = | I - VT VT**T | / ( N ulp ) ', / '11 = 0 if S contains min(M,N) nonnegative values in', ' decreasing order, else 1/ulp', / '12 = | U - Upartial | / ( M ulp )', / '13 = | VT - VTpartial | / ( N ulp )', / '14 = | S - Spartial | / ( min(M,N) ulp |S| )', / '15 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ', / '16 = | I - U**T U | / ( M ulp ) ', / '17 = | I - VT VT**T | / ( N ulp ) ', / '18 = 0 if S contains min(M,N) nonnegative values in', ' decreasing order, else 1/ulp', / '19 = | U - Upartial | / ( M ulp )', / '20 = | VT - VTpartial | / ( N ulp )',/ '21 = | S - Spartial | / ( min(M,N) ulp |S| )',/ '22 = 0 if S contains min(M,N) nonnegative values in',' decreasing order, else 1/ulp',/ '23 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ),',' DGESVDX(V,V,A) ',/ '24 = | I - U**T U | / ( M ulp ) ',/ '25 = | I - VT VT**T | / ( N ulp ) ',/ '26 = 0 if S contains min(M,N) nonnegative values in',' decreasing order, else 1/ulp',/ '27 = | U - Upartial | / ( M ulp )',/ '28 = | VT - VTpartial | / ( N ulp )',/ '29 = | S - Spartial | / ( min(M,N) ulp |S| )',/ '30 = | U**T A VT**T - diag(S) | / ( |A| max(M,N) ulp ),',' DGESVDX(V,V,I) ',/ '31 = | I - U**T U | / ( M ulp ) ',/ '32 = | I - VT VT**T | / ( N ulp ) ',/ '33 = | U**T A VT**T - diag(S) | / ( |A| max(M,N) ulp ),',' DGESVDX(V,V,V) ',/ '34 = | I - U**T U | / ( M ulp ) ',/ '35 = | I - VT VT**T | / ( N ulp ) ',' DGESVDQ(H,N,N,A,A',/ '36 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ',/ '37 = | I - U**T U | / ( M ulp ) ',/ '38 = | I - VT VT**T | / ( N ulp ) ',/ '39 = 0 if S contains min(M,N) nonnegative values in',' decreasing order, else 1/ulp',/ / )
 9997 FORMAT( ' M=', I5, ', N=', I5, ', type ', I1, ', IWS=', I1, ', seed=', 4( I4, ',' ), ' test(', I2, ')=', G11.4 )
 9996 FORMAT( ' DDRVBD: ', A, ' returned INFO=', I6, '.', / 9X, 'M=', I6, ', N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )
 9995 FORMAT( ' DDRVBD: ', A, ' returned INFO=', I6, '.', / 9X, 'M=', I6, ', N=', I6, ', JTYPE=', I6, ', LSWORK=', I6, / 9X, 'ISEED=(', 3( I5, ',' ), I5, ')' )

      RETURN

      // End of DDRVBD

      }
