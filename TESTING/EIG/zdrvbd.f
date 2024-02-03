      SUBROUTINE ZDRVBD( NSIZES, MM, NN, NTYPES, DOTYPE, ISEED, THRESH, A, LDA, U, LDU, VT, LDVT, ASAV, USAV, VTSAV, S, SSAV, E, WORK, LWORK, RWORK, IWORK, NOUNIT, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      IMPLICIT NONE

      // .. Scalar Arguments ..
      int                INFO, LDA, LDU, LDVT, LWORK, NOUNIT, NSIZES, NTYPES;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                ISEED( 4 ), IWORK( * ), MM( * ), NN( * );
      double             E( * ), RWORK( * ), S( * ), SSAV( * );
      COMPLEX*16         A( LDA, * ), ASAV( LDA, * ), U( LDU, * ), USAV( LDU, * ), VT( LDVT, * ), VTSAV( LDVT, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TWO, HALF;
      const              ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, HALF = 0.5D0 ;
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) ;
      int                MAXTYP;
      const              MAXTYP = 5 ;
      // ..
      // .. Local Scalars ..
      bool               BADMM, BADNN;
      String             JOBQ, JOBU, JOBVT, RANGE;
      int                I, IINFO, IJQ, IJU, IJVT, IL, IU, ITEMP, IWSPC, IWTMP, J, JSIZE, JTYPE, LSWORK, M, MINWRK, MMAX, MNMAX, MNMIN, MTYPES, N, NERRS, NFAIL, NMAX, NS, NSI, NSV, NTEST, NTESTF, NTESTT, LRWORK;
      double             ANORM, DIF, DIV, OVFL, RTUNFL, ULP, ULPINV, UNFL, VL, VU;
      // ..
      // .. Local Scalars for ZGESVDQ ..
      int                LIWORK, NUMRANK;
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
      // EXTERNAL ALASVM, XERBLA, ZBDT01, ZBDT05, ZGESDD, ZGESVD, ZGESVDQ, ZGESVJ, ZGEJSV, ZGESVDX, ZLACPY, ZLASET, ZLATMS, ZUNT01, ZUNT03
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, MIN
      // ..
      // .. Scalars in Common ..
      String             SRNAMT;
      // ..
      // .. Common blocks ..
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

      // Important constants

      NERRS = 0
      NTESTT = 0
      NTESTF = 0
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
         MINWRK = MAX( MINWRK, MAX( 3*MIN( MM( J ), NN( J ) )+MAX( MM( J ), NN( J ) )**2, 5*MIN( MM( J ), NN( J ) ), 3*MAX( MM( J ), NN( J ) ) ) )
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
         xerbla('ZDRVBD', -INFO );
         RETURN
      }

      // Quick return if nothing to do

      if (NSIZES == 0 || NTYPES == 0) RETURN;

      // More Important constants

      UNFL = DLAMCH( 'S' )
      OVFL = ONE / UNFL
      ULP = DLAMCH( 'E' )
      ULPINV = ONE / ULP
      RTUNFL = SQRT( UNFL )

      // Loop over sizes, types

      NERRS = 0

      for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) { // 230
         M = MM( JSIZE )
         N = NN( JSIZE )
         MNMIN = MIN( M, N )

         if ( NSIZES != 1 ) {
            MTYPES = MIN( MAXTYP, NTYPES )
         } else {
            MTYPES = MIN( MAXTYP+1, NTYPES )
         }

         for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) { // 220
            IF( .NOT.DOTYPE( JTYPE ) ) GO TO 220
            NTEST = 0

            for (J = 1; J <= 4; J++) { // 20
               IOLDSD( J ) = ISEED( J )
            } // 20

            // Compute "A"

            if (MTYPES > MAXTYP) GO TO 50;

            if ( JTYPE == 1 ) {

               // Zero matrix

               zlaset('Full', M, N, CZERO, CZERO, A, LDA );
               DO 30 I = 1, MIN( M, N )
                  S( I ) = ZERO
               } // 30

            } else if ( JTYPE == 2 ) {

               // Identity matrix

               zlaset('Full', M, N, CZERO, CONE, A, LDA );
               DO 40 I = 1, MIN( M, N )
                  S( I ) = ONE
               } // 40

            } else {

               // (Scaled) random matrix

               if (JTYPE == 3) ANORM = ONE                IF( JTYPE == 4 ) ANORM = UNFL / ULP                IF( JTYPE == 5 ) ANORM = OVFL*ULP;
               zlatms(M, N, 'U', ISEED, 'N', S, 4, DBLE( MNMIN ), ANORM, M-1, N-1, 'N', A, LDA, WORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9996 )'Generator', IINFO, M, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               }
            }

            } // 50
            zlacpy('F', M, N, A, LDA, ASAV, LDA );

            // Do for minimal and adequate (for blocking) workspace

            for (IWSPC = 1; IWSPC <= 4; IWSPC++) { // 210

               // Test for ZGESVD

               IWTMP = 2*MIN( M, N )+MAX( M, N )
               LSWORK = IWTMP + ( IWSPC-1 )*( LWORK-IWTMP ) / 3
               LSWORK = MIN( LSWORK, LWORK )
               LSWORK = MAX( LSWORK, 1 )
               if (IWSPC == 4) LSWORK = LWORK;

               for (J = 1; J <= 35; J++) { // 60
                  RESULT( J ) = -ONE
               } // 60

               // Factorize A

               if (IWSPC > 1) CALL ZLACPY( 'F', M, N, ASAV, LDA, A, LDA );
               SRNAMT = 'ZGESVD'
               zgesvd('A', 'A', M, N, A, LDA, SSAV, USAV, LDU, VTSAV, LDVT, WORK, LSWORK, RWORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9995 )'GESVD', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               }

               // Do tests 1--4

               zbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RWORK, RESULT( 1 ) );
               if ( M != 0 && N != 0 ) {
                  zunt01('Columns', MNMIN, M, USAV, LDU, WORK, LWORK, RWORK, RESULT( 2 ) );
                  zunt01('Rows', MNMIN, N, VTSAV, LDVT, WORK, LWORK, RWORK, RESULT( 3 ) );
               }
               RESULT( 4 ) = 0
               for (I = 1; I <= MNMIN - 1; I++) { // 70
                  IF( SSAV( I ) < SSAV( I+1 ) ) RESULT( 4 ) = ULPINV                   IF( SSAV( I ) < ZERO ) RESULT( 4 ) = ULPINV
               } // 70
               if ( MNMIN >= 1 ) {
                  IF( SSAV( MNMIN ) < ZERO ) RESULT( 4 ) = ULPINV
               }

               // Do partial SVDs, comparing to SSAV, USAV, and VTSAV

               RESULT( 5 ) = ZERO
               RESULT( 6 ) = ZERO
               RESULT( 7 ) = ZERO
               for (IJU = 0; IJU <= 3; IJU++) { // 100
                  for (IJVT = 0; IJVT <= 3; IJVT++) { // 90
                     IF( ( IJU == 3 && IJVT == 3 ) || ( IJU == 1 && IJVT == 1 ) )GO TO 90
                     JOBU = CJOB( IJU+1 )
                     JOBVT = CJOB( IJVT+1 )
                     zlacpy('F', M, N, ASAV, LDA, A, LDA );
                     SRNAMT = 'ZGESVD'
                     zgesvd(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LSWORK, RWORK, IINFO );

                     // Compare U

                     DIF = ZERO
                     if ( M > 0 && N > 0 ) {
                        if ( IJU == 1 ) {
                           zunt03('C', M, MNMIN, M, MNMIN, USAV, LDU, A, LDA, WORK, LWORK, RWORK, DIF, IINFO );
                        } else if ( IJU == 2 ) {
                           zunt03('C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, RWORK, DIF, IINFO );
                        } else if ( IJU == 3 ) {
                           zunt03('C', M, M, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, RWORK, DIF, IINFO );
                        }
                     }
                     RESULT( 5 ) = MAX( RESULT( 5 ), DIF )

                     // Compare VT

                     DIF = ZERO
                     if ( M > 0 && N > 0 ) {
                        if ( IJVT == 1 ) {
                           zunt03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, A, LDA, WORK, LWORK, RWORK, DIF, IINFO );
                        } else if ( IJVT == 2 ) {
                           zunt03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, RWORK, DIF, IINFO );
                        } else if ( IJVT == 3 ) {
                           zunt03('R', N, N, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, RWORK, DIF, IINFO );
                        }
                     }
                     RESULT( 6 ) = MAX( RESULT( 6 ), DIF )

                     // Compare S

                     DIF = ZERO
                     DIV = MAX( DBLE( MNMIN )*ULP*S( 1 ), DLAMCH( 'Safe minimum' ) )
                     for (I = 1; I <= MNMIN - 1; I++) { // 80
                        IF( SSAV( I ) < SSAV( I+1 ) ) DIF = ULPINV                         IF( SSAV( I ) < ZERO ) DIF = ULPINV
                        DIF = MAX( DIF, ABS( SSAV( I )-S( I ) ) / DIV )
                     } // 80
                     RESULT( 7 ) = MAX( RESULT( 7 ), DIF )
                  } // 90
               } // 100

               // Test for ZGESDD

               IWTMP = 2*MNMIN*MNMIN + 2*MNMIN + MAX( M, N )
               LSWORK = IWTMP + ( IWSPC-1 )*( LWORK-IWTMP ) / 3
               LSWORK = MIN( LSWORK, LWORK )
               LSWORK = MAX( LSWORK, 1 )
               if (IWSPC == 4) LSWORK = LWORK;

               // Factorize A

               zlacpy('F', M, N, ASAV, LDA, A, LDA );
               SRNAMT = 'ZGESDD'
               zgesdd('A', M, N, A, LDA, SSAV, USAV, LDU, VTSAV, LDVT, WORK, LSWORK, RWORK, IWORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9995 )'GESDD', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               }

               // Do tests 1--4

               zbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RWORK, RESULT( 8 ) );
               if ( M != 0 && N != 0 ) {
                  zunt01('Columns', MNMIN, M, USAV, LDU, WORK, LWORK, RWORK, RESULT( 9 ) );
                  zunt01('Rows', MNMIN, N, VTSAV, LDVT, WORK, LWORK, RWORK, RESULT( 10 ) );
               }
               RESULT( 11 ) = 0
               for (I = 1; I <= MNMIN - 1; I++) { // 110
                  IF( SSAV( I ) < SSAV( I+1 ) ) RESULT( 11 ) = ULPINV                   IF( SSAV( I ) < ZERO ) RESULT( 11 ) = ULPINV
               } // 110
               if ( MNMIN >= 1 ) {
                  IF( SSAV( MNMIN ) < ZERO ) RESULT( 11 ) = ULPINV
               }

               // Do partial SVDs, comparing to SSAV, USAV, and VTSAV

               RESULT( 12 ) = ZERO
               RESULT( 13 ) = ZERO
               RESULT( 14 ) = ZERO
               for (IJQ = 0; IJQ <= 2; IJQ++) { // 130
                  JOBQ = CJOB( IJQ+1 )
                  zlacpy('F', M, N, ASAV, LDA, A, LDA );
                  SRNAMT = 'ZGESDD'
                  zgesdd(JOBQ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LSWORK, RWORK, IWORK, IINFO );

                  // Compare U

                  DIF = ZERO
                  if ( M > 0 && N > 0 ) {
                     if ( IJQ == 1 ) {
                        if ( M >= N ) {
                           zunt03('C', M, MNMIN, M, MNMIN, USAV, LDU, A, LDA, WORK, LWORK, RWORK, DIF, IINFO );
                        } else {
                           zunt03('C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, RWORK, DIF, IINFO );
                        }
                     } else if ( IJQ == 2 ) {
                        zunt03('C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, RWORK, DIF, IINFO );
                     }
                  }
                  RESULT( 12 ) = MAX( RESULT( 12 ), DIF )

                  // Compare VT

                  DIF = ZERO
                  if ( M > 0 && N > 0 ) {
                     if ( IJQ == 1 ) {
                        if ( M >= N ) {
                           zunt03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, RWORK, DIF, IINFO );
                        } else {
                           zunt03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, A, LDA, WORK, LWORK, RWORK, DIF, IINFO );
                        }
                     } else if ( IJQ == 2 ) {
                        zunt03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, RWORK, DIF, IINFO );
                     }
                  }
                  RESULT( 13 ) = MAX( RESULT( 13 ), DIF )

                  // Compare S

                  DIF = ZERO
                  DIV = MAX( DBLE( MNMIN )*ULP*S( 1 ), DLAMCH( 'Safe minimum' ) )
                  for (I = 1; I <= MNMIN - 1; I++) { // 120
                     IF( SSAV( I ) < SSAV( I+1 ) ) DIF = ULPINV                      IF( SSAV( I ) < ZERO ) DIF = ULPINV
                     DIF = MAX( DIF, ABS( SSAV( I )-S( I ) ) / DIV )
                  } // 120
                  RESULT( 14 ) = MAX( RESULT( 14 ), DIF )
               } // 130

               // Test ZGESVDQ
               // Note: ZGESVDQ only works for M >= N

               RESULT( 36 ) = ZERO
               RESULT( 37 ) = ZERO
               RESULT( 38 ) = ZERO
               RESULT( 39 ) = ZERO

               if ( M >= N ) {
                  IWTMP = 2*MNMIN*MNMIN + 2*MNMIN + MAX( M, N )
                  LSWORK = IWTMP + ( IWSPC-1 )*( LWORK-IWTMP ) / 3
                  LSWORK = MIN( LSWORK, LWORK )
                  LSWORK = MAX( LSWORK, 1 )
                  if (IWSPC == 4) LSWORK = LWORK;

                  zlacpy('F', M, N, ASAV, LDA, A, LDA );
                  SRNAMT = 'ZGESVDQ'

                  LRWORK = MAX(2, M, 5*N)
                  LIWORK = MAX( N, 1 )
                  zgesvdq('H', 'N', 'N', 'A', 'A',  M, N, A, LDA, SSAV, USAV, LDU, VTSAV, LDVT, NUMRANK, IWORK, LIWORK, WORK, LWORK, RWORK, LRWORK, IINFO );

                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9995 )'ZGESVDQ', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                     INFO = ABS( IINFO )
                     RETURN
                  }

                  // Do tests 36--39

                  zbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RWORK, RESULT( 36 ) );
                  if ( M != 0 && N != 0 ) {
                     zunt01('Columns', M, M, USAV, LDU, WORK, LWORK, RWORK, RESULT( 37 ) );
                     zunt01('Rows', N, N, VTSAV, LDVT, WORK, LWORK, RWORK, RESULT( 38 ) );
                  }
                  RESULT( 39 ) = ZERO
                  for (I = 1; I <= MNMIN - 1; I++) { // 199
                     IF( SSAV( I ) < SSAV( I+1 ) ) RESULT( 39 ) = ULPINV                      IF( SSAV( I ) < ZERO ) RESULT( 39 ) = ULPINV
                  } // 199
                  if ( MNMIN >= 1 ) {
                     IF( SSAV( MNMIN ) < ZERO ) RESULT( 39 ) = ULPINV
                  }
               }

               // Test ZGESVJ
               // Note: ZGESVJ only works for M >= N

               RESULT( 15 ) = ZERO
               RESULT( 16 ) = ZERO
               RESULT( 17 ) = ZERO
               RESULT( 18 ) = ZERO

               if ( M >= N ) {
                  IWTMP = 2*MNMIN*MNMIN + 2*MNMIN + MAX( M, N )
                  LSWORK = IWTMP + ( IWSPC-1 )*( LWORK-IWTMP ) / 3
                  LSWORK = MIN( LSWORK, LWORK )
                  LSWORK = MAX( LSWORK, 1 )
                  LRWORK = MAX(6,N)
                  if (IWSPC == 4) LSWORK = LWORK;

                  zlacpy('F', M, N, ASAV, LDA, USAV, LDA );
                  SRNAMT = 'ZGESVJ'
                  zgesvj('G', 'U', 'V', M, N, USAV, LDA, SSAV, 0, A, LDVT, WORK, LWORK, RWORK, LRWORK, IINFO );

                  // ZGESVJ returns V not VH

                  for (J = 1; J <= N; J++) {
                     for (I = 1; I <= N; I++) {
                        VTSAV(J,I) = CONJG (A(I,J))
                     }
                  }

                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9995 )'GESVJ', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                     INFO = ABS( IINFO )
                     RETURN
                  }

                  // Do tests 15--18

                  zbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RWORK, RESULT( 15 ) );
                  if ( M != 0 && N != 0 ) {
                     zunt01('Columns', M, M, USAV, LDU, WORK, LWORK, RWORK, RESULT( 16 ) );
                     zunt01('Rows', N, N, VTSAV, LDVT, WORK, LWORK, RWORK, RESULT( 17 ) );
                  }
                  RESULT( 18 ) = ZERO
                  for (I = 1; I <= MNMIN - 1; I++) { // 131
                     IF( SSAV( I ) < SSAV( I+1 ) ) RESULT( 18 ) = ULPINV                      IF( SSAV( I ) < ZERO ) RESULT( 18 ) = ULPINV
                  } // 131
                  if ( MNMIN >= 1 ) {
                     IF( SSAV( MNMIN ) < ZERO ) RESULT( 18 ) = ULPINV
                  }
               }

               // Test ZGEJSV
               // Note: ZGEJSV only works for M >= N

               RESULT( 19 ) = ZERO
               RESULT( 20 ) = ZERO
               RESULT( 21 ) = ZERO
               RESULT( 22 ) = ZERO
               if ( M >= N ) {
                  IWTMP = 2*MNMIN*MNMIN + 2*MNMIN + MAX( M, N )
                  LSWORK = IWTMP + ( IWSPC-1 )*( LWORK-IWTMP ) / 3
                  LSWORK = MIN( LSWORK, LWORK )
                  LSWORK = MAX( LSWORK, 1 )
                  if (IWSPC == 4) LSWORK = LWORK;
                  LRWORK = MAX( 7, N + 2*M)

                  zlacpy('F', M, N, ASAV, LDA, VTSAV, LDA );
                  SRNAMT = 'ZGEJSV'
                  zgejsv('G', 'U', 'V', 'R', 'N', 'N', M, N, VTSAV, LDA, SSAV, USAV, LDU, A, LDVT, WORK, LWORK, RWORK, LRWORK, IWORK, IINFO );

                  // ZGEJSV returns V not VH

                  for (J = 1; J <= N; J++) { // 133
                     for (I = 1; I <= N; I++) { // 132
                        VTSAV(J,I) = CONJG (A(I,J))
  132                END DO
  133             END DO

                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9995 )'GEJSV', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                     INFO = ABS( IINFO )
                     RETURN
                  }

                  // Do tests 19--22

                  zbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RWORK, RESULT( 19 ) );
                  if ( M != 0 && N != 0 ) {
                     zunt01('Columns', M, M, USAV, LDU, WORK, LWORK, RWORK, RESULT( 20 ) );
                     zunt01('Rows', N, N, VTSAV, LDVT, WORK, LWORK, RWORK, RESULT( 21 ) );
                  }
                  RESULT( 22 ) = ZERO
                  for (I = 1; I <= MNMIN - 1; I++) { // 134
                     IF( SSAV( I ) < SSAV( I+1 ) ) RESULT( 22 ) = ULPINV                      IF( SSAV( I ) < ZERO ) RESULT( 22 ) = ULPINV
                  } // 134
                  if ( MNMIN >= 1 ) {
                     IF( SSAV( MNMIN ) < ZERO ) RESULT( 22 ) = ULPINV
                  }
               }

               // Test ZGESVDX

               // Factorize A

               zlacpy('F', M, N, ASAV, LDA, A, LDA );
               SRNAMT = 'ZGESVDX'
               zgesvdx('V', 'V', 'A', M, N, A, LDA, VL, VU, IL, IU, NS, SSAV, USAV, LDU, VTSAV, LDVT, WORK, LWORK, RWORK, IWORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9995 )'GESVDX', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               }

               // Do tests 1--4

               RESULT( 23 ) = ZERO
               RESULT( 24 ) = ZERO
               RESULT( 25 ) = ZERO
               zbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RWORK, RESULT( 23 ) );
               if ( M != 0 && N != 0 ) {
                  zunt01('Columns', MNMIN, M, USAV, LDU, WORK, LWORK, RWORK, RESULT( 24 ) );
                  zunt01('Rows', MNMIN, N, VTSAV, LDVT, WORK, LWORK, RWORK, RESULT( 25 ) );
               }
               RESULT( 26 ) = ZERO
               for (I = 1; I <= MNMIN - 1; I++) { // 140
                  IF( SSAV( I ) < SSAV( I+1 ) ) RESULT( 26 ) = ULPINV                   IF( SSAV( I ) < ZERO ) RESULT( 26 ) = ULPINV
               } // 140
               if ( MNMIN >= 1 ) {
                  IF( SSAV( MNMIN ) < ZERO ) RESULT( 26 ) = ULPINV
               }

               // Do partial SVDs, comparing to SSAV, USAV, and VTSAV

               RESULT( 27 ) = ZERO
               RESULT( 28 ) = ZERO
               RESULT( 29 ) = ZERO
               for (IJU = 0; IJU <= 1; IJU++) { // 170
                  for (IJVT = 0; IJVT <= 1; IJVT++) { // 160
                     IF( ( IJU == 0 && IJVT == 0 ) || ( IJU == 1 && IJVT == 1 ) ) GO TO 160
                     JOBU = CJOBV( IJU+1 )
                     JOBVT = CJOBV( IJVT+1 )
                     RANGE = CJOBR( 1 )
                     zlacpy('F', M, N, ASAV, LDA, A, LDA );
                     SRNAMT = 'ZGESVDX'
                     zgesvdx(JOBU, JOBVT, 'A', M, N, A, LDA, VL, VU, IL, IU, NS, SSAV, U, LDU, VT, LDVT, WORK, LWORK, RWORK, IWORK, IINFO );

                     // Compare U

                     DIF = ZERO
                     if ( M > 0 && N > 0 ) {
                        if ( IJU == 1 ) {
                           zunt03('C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, RWORK, DIF, IINFO );
                        }
                     }
                     RESULT( 27 ) = MAX( RESULT( 27 ), DIF )

                     // Compare VT

                     DIF = ZERO
                     if ( M > 0 && N > 0 ) {
                        if ( IJVT == 1 ) {
                           zunt03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, RWORK, DIF, IINFO );
                        }
                     }
                     RESULT( 28 ) = MAX( RESULT( 28 ), DIF )

                     // Compare S

                     DIF = ZERO
                     DIV = MAX( DBLE( MNMIN )*ULP*S( 1 ), DLAMCH( 'Safe minimum' ) )
                     for (I = 1; I <= MNMIN - 1; I++) { // 150
                        IF( SSAV( I ) < SSAV( I+1 ) ) DIF = ULPINV                         IF( SSAV( I ) < ZERO ) DIF = ULPINV
                        DIF = MAX( DIF, ABS( SSAV( I )-S( I ) ) / DIV )
                     } // 150
                     RESULT( 29) = MAX( RESULT( 29 ), DIF )
                  } // 160
               } // 170

               // Do tests 8--10

               for (I = 1; I <= 4; I++) { // 180
                  ISEED2( I ) = ISEED( I )
               } // 180
               if ( MNMIN.LE.1 ) {
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
               zlacpy('F', M, N, ASAV, LDA, A, LDA );
               SRNAMT = 'ZGESVDX'
               zgesvdx('V', 'V', 'I', M, N, A, LDA, VL, VU, IL, IU, NSI, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, IWORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9995 )'GESVDX', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               }

               RESULT( 30 ) = ZERO
               RESULT( 31 ) = ZERO
               RESULT( 32 ) = ZERO
               zbdt05(M, N, ASAV, LDA, S, NSI, U, LDU, VT, LDVT, WORK, RESULT( 30 ) );
               if ( M != 0 && N != 0 ) {
                  zunt01('Columns', M, NSI, U, LDU, WORK, LWORK, RWORK, RESULT( 31 ) );
                  zunt01('Rows', NSI, N, VT, LDVT, WORK, LWORK, RWORK, RESULT( 32 ) );
               }

               // Do tests 11--13

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
               zlacpy('F', M, N, ASAV, LDA, A, LDA );
               SRNAMT = 'ZGESVDX'
               zgesvdx('V', 'V', 'V', M, N, A, LDA, VL, VU, IL, IU, NSV, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, IWORK, IINFO );
               if ( IINFO != 0 ) {
                  WRITE( NOUNIT, FMT = 9995 )'GESVDX', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               }

               RESULT( 33 ) = ZERO
               RESULT( 34 ) = ZERO
               RESULT( 35 ) = ZERO
               zbdt05(M, N, ASAV, LDA, S, NSV, U, LDU, VT, LDVT, WORK, RESULT( 33 ) );
               if ( M != 0 && N != 0 ) {
                  zunt01('Columns', M, NSV, U, LDU, WORK, LWORK, RWORK, RESULT( 34 ) );
                  zunt01('Rows', NSV, N, VT, LDVT, WORK, LWORK, RWORK, RESULT( 35 ) );
               }

               // End of Loop -- Check for RESULT(j) > THRESH

               NTEST = 0
               NFAIL = 0
               for (J = 1; J <= 39; J++) { // 190
                  IF( RESULT( J ) >= ZERO ) NTEST = NTEST + 1                   IF( RESULT( J ) >= THRESH ) NFAIL = NFAIL + 1
               } // 190

               if (NFAIL > 0) NTESTF = NTESTF + 1;
               if ( NTESTF == 1 ) {
                  WRITE( NOUNIT, FMT = 9999 )
                  WRITE( NOUNIT, FMT = 9998 )THRESH
                  NTESTF = 2
               }

               for (J = 1; J <= 39; J++) { // 200
                  if ( RESULT( J ) >= THRESH ) {
                     WRITE( NOUNIT, FMT = 9997 )M, N, JTYPE, IWSPC, IOLDSD, J, RESULT( J )
                  }
               } // 200

               NERRS = NERRS + NFAIL
               NTESTT = NTESTT + NTEST

            } // 210

         } // 220
      } // 230

      // Summary

      alasvm('ZBD', NOUNIT, NERRS, NTESTT, 0 );

 9999 FORMAT( ' SVD -- Complex Singular Value Decomposition Driver ', / ' Matrix types (see ZDRVBD for details):', / / ' 1 = Zero matrix', / ' 2 = Identity matrix', / ' 3 = Evenly spaced singular values near 1', / ' 4 = Evenly spaced singular values near underflow', / ' 5 = Evenly spaced singular values near overflow', / / ' Tests performed: ( A is dense, U and V are unitary,', / 19X, ' S is an array, and Upartial, VTpartial, and', / 19X, ' Spartial are partially computed U, VT and S),', / )
 9998 FORMAT( ' Tests performed with Test Threshold = ', F8.2, / ' ZGESVD: ', / ' 1 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ', / ' 2 = | I - U**T U | / ( M ulp ) ', / ' 3 = | I - VT VT**T | / ( N ulp ) ', / ' 4 = 0 if S contains min(M,N) nonnegative values in', ' decreasing order, else 1/ulp', / ' 5 = | U - Upartial | / ( M ulp )', / ' 6 = | VT - VTpartial | / ( N ulp )', / ' 7 = | S - Spartial | / ( min(M,N) ulp |S| )', / ' ZGESDD: ', / ' 8 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ', / ' 9 = | I - U**T U | / ( M ulp ) ', / '10 = | I - VT VT**T | / ( N ulp ) ', / '11 = 0 if S contains min(M,N) nonnegative values in', ' decreasing order, else 1/ulp', / '12 = | U - Upartial | / ( M ulp )', / '13 = | VT - VTpartial | / ( N ulp )', / '14 = | S - Spartial | / ( min(M,N) ulp |S| )', / ' ZGESVJ: ', / / '15 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ', / '16 = | I - U**T U | / ( M ulp ) ', / '17 = | I - VT VT**T | / ( N ulp ) ',/ '18 = 0 if S contains min(M,N) nonnegative values in',' decreasing order, else 1/ulp',/ ' ZGESJV: ', // '19 = | A - U diag(S) VT | / ( |A| max(M,N) ulp )',/ '20 = | I - U**T U | / ( M ulp ) ',/ '21 = | I - VT VT**T | / ( N ulp ) ',/ '22 = 0 if S contains min(M,N) nonnegative values in',' decreasing order, else 1/ulp',/ ' ZGESVDX(V,V,A): ', /  '23 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ',/ '24 = | I - U**T U | / ( M ulp ) ',/ '25 = | I - VT VT**T | / ( N ulp ) ',/ '26 = 0 if S contains min(M,N) nonnegative values in',' decreasing order, else 1/ulp',/ '27 = | U - Upartial | / ( M ulp )',/ '28 = | VT - VTpartial | / ( N ulp )',/ '29 = | S - Spartial | / ( min(M,N) ulp |S| )',/ ' ZGESVDX(V,V,I): ',/ '30 = | U**T A VT**T - diag(S) | / ( |A| max(M,N) ulp )',/ '31 = | I - U**T U | / ( M ulp ) ',/ '32 = | I - VT VT**T | / ( N ulp ) ',/ ' ZGESVDX(V,V,V) ',/ '33 = | U**T A VT**T - diag(S) | / ( |A| max(M,N) ulp )',/ '34 = | I - U**T U | / ( M ulp ) ',/ '35 = | I - VT VT**T | / ( N ulp ) ',' ZGESVDQ(H,N,N,A,A',/ '36 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ',/ '37 = | I - U**T U | / ( M ulp ) ',/ '38 = | I - VT VT**T | / ( N ulp ) ',/ '39 = 0 if S contains min(M,N) nonnegative values in',' decreasing order, else 1/ulp',/ / )
 9997 FORMAT( ' M=', I5, ', N=', I5, ', type ', I1, ', IWS=', I1, ', seed=', 4( I4, ',' ), ' test(', I2, ')=', G11.4 )
 9996 FORMAT( ' ZDRVBD: ', A, ' returned INFO=', I6, '.', / 9X, 'M=', I6, ', N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )
 9995 FORMAT( ' ZDRVBD: ', A, ' returned INFO=', I6, '.', / 9X, 'M=', I6, ', N=', I6, ', JTYPE=', I6, ', LSWORK=', I6, / 9X, 'ISEED=(', 3( I5, ',' ), I5, ')' )

      RETURN

      // End of ZDRVBD

      }
