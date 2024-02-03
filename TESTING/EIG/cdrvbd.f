      SUBROUTINE CDRVBD( NSIZES, MM, NN, NTYPES, DOTYPE, ISEED, THRESH, A, LDA, U, LDU, VT, LDVT, ASAV, USAV, VTSAV, S, SSAV, E, WORK, LWORK, RWORK, IWORK, NOUNIT, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      IMPLICIT NONE

      // .. Scalar Arguments ..
      int                INFO, LDA, LDU, LDVT, LWORK, NOUNIT, NSIZES, NTYPES;
      REAL               THRESH
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                ISEED( 4 ), IWORK( * ), MM( * ), NN( * );
      REAL               E( * ), RWORK( * ), S( * ), SSAV( * )
      COMPLEX            A( LDA, * ), ASAV( LDA, * ), U( LDU, * ), USAV( LDU, * ), VT( LDVT, * ), VTSAV( LDVT, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TWO, HALF
      const              ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0, HALF = 0.5E0 ;
      COMPLEX            CZERO, CONE
      const              CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) ;
      int                MAXTYP;
      const              MAXTYP = 5 ;
      // ..
      // .. Local Scalars ..
      bool               BADMM, BADNN;
      String             JOBQ, JOBU, JOBVT, RANGE;
      int                I, IINFO, IJQ, IJU, IJVT, IL, IU, ITEMP, IWSPC, IWTMP, J, JSIZE, JTYPE, LSWORK, M, MINWRK, MMAX, MNMAX, MNMIN, MTYPES, N, NERRS, NFAIL, NMAX, NS, NSI, NSV, NTEST, NTESTF, NTESTT, LRWORK;
      REAL               ANORM, DIF, DIV, OVFL, RTUNFL, ULP, ULPINV, UNFL, VL, VU
      // ..
      // .. Local Scalars for CGESVDQ ..
      int                LIWORK, NUMRANK;
      // ..
      // .. Local Arrays ..
      String             CJOB( 4 ), CJOBR( 3 ), CJOBV( 2 );
      int                IOLDSD( 4 ), ISEED2( 4 );
      REAL               RESULT( 39 )
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SLARND
      // EXTERNAL SLAMCH, SLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALASVM, XERBLA, CBDT01, CBDT05, CGESDD, CGESVD, CGESVDQ, CGESVJ, CGEJSV, CGESVDX, CLACPY, CLASET, CLATMS, CUNT01, CUNT03
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, REAL, MAX, MIN
      // ..
      // .. Scalars in Common ..
      String             SRNAMT;
      // ..
      // .. Common blocks ..
      COMMON             / SRNAMC / SRNAMT
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
      BADMM = .FALSE.
      BADNN = .FALSE.
      MMAX = 1
      NMAX = 1
      MNMAX = 1
      MINWRK = 1
      DO 10 J = 1, NSIZES
         MMAX = MAX( MMAX, MM( J ) )
         IF( MM( J ).LT.0 ) BADMM = .TRUE.
         NMAX = MAX( NMAX, NN( J ) )
         IF( NN( J ).LT.0 ) BADNN = .TRUE.
         MNMAX = MAX( MNMAX, MIN( MM( J ), NN( J ) ) )
         MINWRK = MAX( MINWRK, MAX( 3*MIN( MM( J ), NN( J ) )+MAX( MM( J ), NN( J ) )**2, 5*MIN( MM( J ), NN( J ) ), 3*MAX( MM( J ), NN( J ) ) ) )
   10 CONTINUE

      // Check for errors

      if ( NSIZES.LT.0 ) {
         INFO = -1
      } else if ( BADMM ) {
         INFO = -2
      } else if ( BADNN ) {
         INFO = -3
      } else if ( NTYPES.LT.0 ) {
         INFO = -4
      } else if ( LDA.LT.MAX( 1, MMAX ) ) {
         INFO = -10
      } else if ( LDU.LT.MAX( 1, MMAX ) ) {
         INFO = -12
      } else if ( LDVT.LT.MAX( 1, NMAX ) ) {
         INFO = -14
      } else if ( MINWRK.GT.LWORK ) {
         INFO = -21
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'CDRVBD', -INFO )
         RETURN
      }

      // Quick return if nothing to do

      IF( NSIZES.EQ.0 .OR. NTYPES.EQ.0 ) RETURN

      // More Important constants

      UNFL = SLAMCH( 'S' )
      OVFL = ONE / UNFL
      ULP = SLAMCH( 'E' )
      ULPINV = ONE / ULP
      RTUNFL = SQRT( UNFL )

      // Loop over sizes, types

      NERRS = 0

      DO 310 JSIZE = 1, NSIZES
         M = MM( JSIZE )
         N = NN( JSIZE )
         MNMIN = MIN( M, N )

         if ( NSIZES.NE.1 ) {
            MTYPES = MIN( MAXTYP, NTYPES )
         } else {
            MTYPES = MIN( MAXTYP+1, NTYPES )
         }

         DO 300 JTYPE = 1, MTYPES
            IF( .NOT.DOTYPE( JTYPE ) ) GO TO 300
            NTEST = 0

            DO 20 J = 1, 4
               IOLDSD( J ) = ISEED( J )
   20       CONTINUE

            // Compute "A"

            IF( MTYPES.GT.MAXTYP ) GO TO 50

            if ( JTYPE.EQ.1 ) {

               // Zero matrix

               CALL CLASET( 'Full', M, N, CZERO, CZERO, A, LDA )
               DO 30 I = 1, MIN( M, N )
                  S( I ) = ZERO
   30          CONTINUE

            } else if ( JTYPE.EQ.2 ) {

               // Identity matrix

               CALL CLASET( 'Full', M, N, CZERO, CONE, A, LDA )
               DO 40 I = 1, MIN( M, N )
                  S( I ) = ONE
   40          CONTINUE

            } else {

               // (Scaled) random matrix

               IF( JTYPE.EQ.3 ) ANORM = ONE                IF( JTYPE.EQ.4 ) ANORM = UNFL / ULP                IF( JTYPE.EQ.5 ) ANORM = OVFL*ULP                CALL CLATMS( M, N, 'U', ISEED, 'N', S, 4, REAL( MNMIN ), ANORM, M-1, N-1, 'N', A, LDA, WORK, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9996 )'Generator', IINFO, M, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               }
            }

   50       CONTINUE
            CALL CLACPY( 'F', M, N, A, LDA, ASAV, LDA )

            // Do for minimal and adequate (for blocking) workspace

            DO 290 IWSPC = 1, 4

               // Test for CGESVD

               IWTMP = 2*MIN( M, N )+MAX( M, N )
               LSWORK = IWTMP + ( IWSPC-1 )*( LWORK-IWTMP ) / 3
               LSWORK = MIN( LSWORK, LWORK )
               LSWORK = MAX( LSWORK, 1 )
               IF( IWSPC.EQ.4 ) LSWORK = LWORK

               DO 60 J = 1, 35
                  RESULT( J ) = -ONE
   60          CONTINUE

               // Factorize A

               IF( IWSPC.GT.1 ) CALL CLACPY( 'F', M, N, ASAV, LDA, A, LDA )
               SRNAMT = 'CGESVD'
               CALL CGESVD( 'A', 'A', M, N, A, LDA, SSAV, USAV, LDU, VTSAV, LDVT, WORK, LSWORK, RWORK, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9995 )'GESVD', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               }

               // Do tests 1--4

               CALL CBDT01( M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RWORK, RESULT( 1 ) )
               if ( M.NE.0 .AND. N.NE.0 ) {
                  CALL CUNT01( 'Columns', MNMIN, M, USAV, LDU, WORK, LWORK, RWORK, RESULT( 2 ) )                   CALL CUNT01( 'Rows', MNMIN, N, VTSAV, LDVT, WORK, LWORK, RWORK, RESULT( 3 ) )
               }
               RESULT( 4 ) = 0
               DO 70 I = 1, MNMIN - 1
                  IF( SSAV( I ).LT.SSAV( I+1 ) ) RESULT( 4 ) = ULPINV                   IF( SSAV( I ).LT.ZERO ) RESULT( 4 ) = ULPINV
   70          CONTINUE
               if ( MNMIN.GE.1 ) {
                  IF( SSAV( MNMIN ).LT.ZERO ) RESULT( 4 ) = ULPINV
               }

               // Do partial SVDs, comparing to SSAV, USAV, and VTSAV

               RESULT( 5 ) = ZERO
               RESULT( 6 ) = ZERO
               RESULT( 7 ) = ZERO
               DO 100 IJU = 0, 3
                  DO 90 IJVT = 0, 3
                     IF( ( IJU.EQ.3 .AND. IJVT.EQ.3 ) .OR. ( IJU.EQ.1 .AND. IJVT.EQ.1 ) )GO TO 90
                     JOBU = CJOB( IJU+1 )
                     JOBVT = CJOB( IJVT+1 )
                     CALL CLACPY( 'F', M, N, ASAV, LDA, A, LDA )
                     SRNAMT = 'CGESVD'
                     CALL CGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LSWORK, RWORK, IINFO )

                     // Compare U

                     DIF = ZERO
                     if ( M.GT.0 .AND. N.GT.0 ) {
                        if ( IJU.EQ.1 ) {
                           CALL CUNT03( 'C', M, MNMIN, M, MNMIN, USAV, LDU, A, LDA, WORK, LWORK, RWORK, DIF, IINFO )
                        } else if ( IJU.EQ.2 ) {
                           CALL CUNT03( 'C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, RWORK, DIF, IINFO )
                        } else if ( IJU.EQ.3 ) {
                           CALL CUNT03( 'C', M, M, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, RWORK, DIF, IINFO )
                        }
                     }
                     RESULT( 5 ) = MAX( RESULT( 5 ), DIF )

                     // Compare VT

                     DIF = ZERO
                     if ( M.GT.0 .AND. N.GT.0 ) {
                        if ( IJVT.EQ.1 ) {
                           CALL CUNT03( 'R', N, MNMIN, N, MNMIN, VTSAV, LDVT, A, LDA, WORK, LWORK, RWORK, DIF, IINFO )
                        } else if ( IJVT.EQ.2 ) {
                           CALL CUNT03( 'R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, RWORK, DIF, IINFO )
                        } else if ( IJVT.EQ.3 ) {
                           CALL CUNT03( 'R', N, N, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, RWORK, DIF, IINFO )
                        }
                     }
                     RESULT( 6 ) = MAX( RESULT( 6 ), DIF )

                     // Compare S

                     DIF = ZERO
                     DIV = MAX( REAL( MNMIN )*ULP*S( 1 ), SLAMCH( 'Safe minimum' ) )
                     DO 80 I = 1, MNMIN - 1
                        IF( SSAV( I ).LT.SSAV( I+1 ) ) DIF = ULPINV                         IF( SSAV( I ).LT.ZERO ) DIF = ULPINV
                        DIF = MAX( DIF, ABS( SSAV( I )-S( I ) ) / DIV )
   80                CONTINUE
                     RESULT( 7 ) = MAX( RESULT( 7 ), DIF )
   90             CONTINUE
  100          CONTINUE

               // Test for CGESDD

               IWTMP = 2*MNMIN*MNMIN + 2*MNMIN + MAX( M, N )
               LSWORK = IWTMP + ( IWSPC-1 )*( LWORK-IWTMP ) / 3
               LSWORK = MIN( LSWORK, LWORK )
               LSWORK = MAX( LSWORK, 1 )
               IF( IWSPC.EQ.4 ) LSWORK = LWORK

               // Factorize A

               CALL CLACPY( 'F', M, N, ASAV, LDA, A, LDA )
               SRNAMT = 'CGESDD'
               CALL CGESDD( 'A', M, N, A, LDA, SSAV, USAV, LDU, VTSAV, LDVT, WORK, LSWORK, RWORK, IWORK, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9995 )'GESDD', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               }

               // Do tests 1--4

               CALL CBDT01( M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RWORK, RESULT( 8 ) )
               if ( M.NE.0 .AND. N.NE.0 ) {
                  CALL CUNT01( 'Columns', MNMIN, M, USAV, LDU, WORK, LWORK, RWORK, RESULT( 9 ) )                   CALL CUNT01( 'Rows', MNMIN, N, VTSAV, LDVT, WORK, LWORK, RWORK, RESULT( 10 ) )
               }
               RESULT( 11 ) = 0
               DO 110 I = 1, MNMIN - 1
                  IF( SSAV( I ).LT.SSAV( I+1 ) ) RESULT( 11 ) = ULPINV                   IF( SSAV( I ).LT.ZERO ) RESULT( 11 ) = ULPINV
  110          CONTINUE
               if ( MNMIN.GE.1 ) {
                  IF( SSAV( MNMIN ).LT.ZERO ) RESULT( 11 ) = ULPINV
               }

               // Do partial SVDs, comparing to SSAV, USAV, and VTSAV

               RESULT( 12 ) = ZERO
               RESULT( 13 ) = ZERO
               RESULT( 14 ) = ZERO
               DO 130 IJQ = 0, 2
                  JOBQ = CJOB( IJQ+1 )
                  CALL CLACPY( 'F', M, N, ASAV, LDA, A, LDA )
                  SRNAMT = 'CGESDD'
                  CALL CGESDD( JOBQ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LSWORK, RWORK, IWORK, IINFO )

                  // Compare U

                  DIF = ZERO
                  if ( M.GT.0 .AND. N.GT.0 ) {
                     if ( IJQ.EQ.1 ) {
                        if ( M.GE.N ) {
                           CALL CUNT03( 'C', M, MNMIN, M, MNMIN, USAV, LDU, A, LDA, WORK, LWORK, RWORK, DIF, IINFO )
                        } else {
                           CALL CUNT03( 'C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, RWORK, DIF, IINFO )
                        }
                     } else if ( IJQ.EQ.2 ) {
                        CALL CUNT03( 'C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, RWORK, DIF, IINFO )
                     }
                  }
                  RESULT( 12 ) = MAX( RESULT( 12 ), DIF )

                  // Compare VT

                  DIF = ZERO
                  if ( M.GT.0 .AND. N.GT.0 ) {
                     if ( IJQ.EQ.1 ) {
                        if ( M.GE.N ) {
                           CALL CUNT03( 'R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, RWORK, DIF, IINFO )
                        } else {
                           CALL CUNT03( 'R', N, MNMIN, N, MNMIN, VTSAV, LDVT, A, LDA, WORK, LWORK, RWORK, DIF, IINFO )
                        }
                     } else if ( IJQ.EQ.2 ) {
                        CALL CUNT03( 'R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, RWORK, DIF, IINFO )
                     }
                  }
                  RESULT( 13 ) = MAX( RESULT( 13 ), DIF )

                  // Compare S

                  DIF = ZERO
                  DIV = MAX( REAL( MNMIN )*ULP*S( 1 ), SLAMCH( 'Safe minimum' ) )
                  DO 120 I = 1, MNMIN - 1
                     IF( SSAV( I ).LT.SSAV( I+1 ) ) DIF = ULPINV                      IF( SSAV( I ).LT.ZERO ) DIF = ULPINV
                     DIF = MAX( DIF, ABS( SSAV( I )-S( I ) ) / DIV )
  120             CONTINUE
                  RESULT( 14 ) = MAX( RESULT( 14 ), DIF )
  130          CONTINUE


               // Test CGESVDQ
               // Note: CGESVDQ only works for M >= N

               RESULT( 36 ) = ZERO
               RESULT( 37 ) = ZERO
               RESULT( 38 ) = ZERO
               RESULT( 39 ) = ZERO

               if ( M.GE.N ) {
                  IWTMP = 2*MNMIN*MNMIN + 2*MNMIN + MAX( M, N )
                  LSWORK = IWTMP + ( IWSPC-1 )*( LWORK-IWTMP ) / 3
                  LSWORK = MIN( LSWORK, LWORK )
                  LSWORK = MAX( LSWORK, 1 )
                  IF( IWSPC.EQ.4 ) LSWORK = LWORK

                  CALL CLACPY( 'F', M, N, ASAV, LDA, A, LDA )
                  SRNAMT = 'CGESVDQ'

                  LRWORK = MAX(2, M, 5*N)
                  LIWORK = MAX( N, 1 )
                  CALL CGESVDQ( 'H', 'N', 'N', 'A', 'A',  M, N, A, LDA, SSAV, USAV, LDU, VTSAV, LDVT, NUMRANK, IWORK, LIWORK, WORK, LWORK, RWORK, LRWORK, IINFO )

                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9995 )'CGESVDQ', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                     INFO = ABS( IINFO )
                     RETURN
                  }

                  // Do tests 36--39

                  CALL CBDT01( M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RWORK, RESULT( 36 ) )
                  if ( M.NE.0 .AND. N.NE.0 ) {
                     CALL CUNT01( 'Columns', M, M, USAV, LDU, WORK, LWORK, RWORK, RESULT( 37 ) )                      CALL CUNT01( 'Rows', N, N, VTSAV, LDVT, WORK, LWORK, RWORK, RESULT( 38 ) )
                  }
                  RESULT( 39 ) = ZERO
                  DO 199 I = 1, MNMIN - 1
                     IF( SSAV( I ).LT.SSAV( I+1 ) ) RESULT( 39 ) = ULPINV                      IF( SSAV( I ).LT.ZERO ) RESULT( 39 ) = ULPINV
  199             CONTINUE
                  if ( MNMIN.GE.1 ) {
                     IF( SSAV( MNMIN ).LT.ZERO ) RESULT( 39 ) = ULPINV
                  }
               }

               // Test CGESVJ
               // Note: CGESVJ only works for M >= N

               RESULT( 15 ) = ZERO
               RESULT( 16 ) = ZERO
               RESULT( 17 ) = ZERO
               RESULT( 18 ) = ZERO

               if ( M.GE.N ) {
                  IWTMP = 2*MNMIN*MNMIN + 2*MNMIN + MAX( M, N )
                  LSWORK = IWTMP + ( IWSPC-1 )*( LWORK-IWTMP ) / 3
                  LSWORK = MIN( LSWORK, LWORK )
                  LSWORK = MAX( LSWORK, 1 )
                  LRWORK = MAX(6,N)
                  IF( IWSPC.EQ.4 ) LSWORK = LWORK

                  CALL CLACPY( 'F', M, N, ASAV, LDA, USAV, LDA )
                  SRNAMT = 'CGESVJ'
                  CALL CGESVJ( 'G', 'U', 'V', M, N, USAV, LDA, SSAV, 0, A, LDVT, WORK, LWORK, RWORK, LRWORK, IINFO )

                  // CGESVJ returns V not VH

                  DO J=1,N
                     DO I=1,N
                        VTSAV(J,I) = CONJG (A(I,J))
                     END DO
                  END DO

                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9995 )'GESVJ', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                     INFO = ABS( IINFO )
                     RETURN
                  }

                  // Do tests 15--18

                  CALL CBDT01( M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RWORK, RESULT( 15 ) )
                  if ( M.NE.0 .AND. N.NE.0 ) {
                     CALL CUNT01( 'Columns', M, M, USAV, LDU, WORK, LWORK, RWORK, RESULT( 16 ) )                      CALL CUNT01( 'Rows', N, N, VTSAV, LDVT, WORK, LWORK, RWORK, RESULT( 17 ) )
                  }
                  RESULT( 18 ) = ZERO
                  DO 131 I = 1, MNMIN - 1
                     IF( SSAV( I ).LT.SSAV( I+1 ) ) RESULT( 18 ) = ULPINV                      IF( SSAV( I ).LT.ZERO ) RESULT( 18 ) = ULPINV
  131             CONTINUE
                  if ( MNMIN.GE.1 ) {
                     IF( SSAV( MNMIN ).LT.ZERO ) RESULT( 18 ) = ULPINV
                  }
               }

               // Test CGEJSV
               // Note: CGEJSV only works for M >= N

               RESULT( 19 ) = ZERO
               RESULT( 20 ) = ZERO
               RESULT( 21 ) = ZERO
               RESULT( 22 ) = ZERO
               if ( M.GE.N ) {
                  IWTMP = 2*MNMIN*MNMIN + 2*MNMIN + MAX( M, N )
                  LSWORK = IWTMP + ( IWSPC-1 )*( LWORK-IWTMP ) / 3
                  LSWORK = MIN( LSWORK, LWORK )
                  LSWORK = MAX( LSWORK, 1 )
                  IF( IWSPC.EQ.4 ) LSWORK = LWORK
                  LRWORK = MAX( 7, N + 2*M)

                  CALL CLACPY( 'F', M, N, ASAV, LDA, VTSAV, LDA )
                  SRNAMT = 'CGEJSV'
                  CALL CGEJSV( 'G', 'U', 'V', 'R', 'N', 'N', M, N, VTSAV, LDA, SSAV, USAV, LDU, A, LDVT, WORK, LWORK, RWORK, LRWORK, IWORK, IINFO )

                  // CGEJSV returns V not VH

                  DO 133 J=1,N
                     DO 132 I=1,N
                        VTSAV(J,I) = CONJG (A(I,J))
  132                END DO
  133             END DO

                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9995 )'GEJSV', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                     INFO = ABS( IINFO )
                     RETURN
                  }

                  // Do tests 19--22

                  CALL CBDT01( M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RWORK, RESULT( 19 ) )
                  if ( M.NE.0 .AND. N.NE.0 ) {
                     CALL CUNT01( 'Columns', M, M, USAV, LDU, WORK, LWORK, RWORK, RESULT( 20 ) )                      CALL CUNT01( 'Rows', N, N, VTSAV, LDVT, WORK, LWORK, RWORK, RESULT( 21 ) )
                  }
                  RESULT( 22 ) = ZERO
                  DO 134 I = 1, MNMIN - 1
                     IF( SSAV( I ).LT.SSAV( I+1 ) ) RESULT( 22 ) = ULPINV                      IF( SSAV( I ).LT.ZERO ) RESULT( 22 ) = ULPINV
  134             CONTINUE
                  if ( MNMIN.GE.1 ) {
                     IF( SSAV( MNMIN ).LT.ZERO ) RESULT( 22 ) = ULPINV
                  }
               }

               // Test CGESVDX

               // Factorize A

               CALL CLACPY( 'F', M, N, ASAV, LDA, A, LDA )
               SRNAMT = 'CGESVDX'
               CALL CGESVDX( 'V', 'V', 'A', M, N, A, LDA, VL, VU, IL, IU, NS, SSAV, USAV, LDU, VTSAV, LDVT, WORK, LWORK, RWORK, IWORK, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9995 )'GESVDX', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               }

               // Do tests 1--4

               RESULT( 23 ) = ZERO
               RESULT( 24 ) = ZERO
               RESULT( 25 ) = ZERO
               CALL CBDT01( M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RWORK, RESULT( 23 ) )
               if ( M.NE.0 .AND. N.NE.0 ) {
                  CALL CUNT01( 'Columns', MNMIN, M, USAV, LDU, WORK, LWORK, RWORK, RESULT( 24 ) )                   CALL CUNT01( 'Rows', MNMIN, N, VTSAV, LDVT, WORK, LWORK, RWORK, RESULT( 25 ) )
               }
               RESULT( 26 ) = ZERO
               DO 140 I = 1, MNMIN - 1
                  IF( SSAV( I ).LT.SSAV( I+1 ) ) RESULT( 26 ) = ULPINV                   IF( SSAV( I ).LT.ZERO ) RESULT( 26 ) = ULPINV
  140          CONTINUE
               if ( MNMIN.GE.1 ) {
                  IF( SSAV( MNMIN ).LT.ZERO ) RESULT( 26 ) = ULPINV
               }

               // Do partial SVDs, comparing to SSAV, USAV, and VTSAV

               RESULT( 27 ) = ZERO
               RESULT( 28 ) = ZERO
               RESULT( 29 ) = ZERO
               DO 170 IJU = 0, 1
                  DO 160 IJVT = 0, 1
                     IF( ( IJU.EQ.0 .AND. IJVT.EQ.0 ) .OR. ( IJU.EQ.1 .AND. IJVT.EQ.1 ) ) GO TO 160
                     JOBU = CJOBV( IJU+1 )
                     JOBVT = CJOBV( IJVT+1 )
                     RANGE = CJOBR( 1 )
                     CALL CLACPY( 'F', M, N, ASAV, LDA, A, LDA )
                     SRNAMT = 'CGESVDX'
                     CALL CGESVDX( JOBU, JOBVT, 'A', M, N, A, LDA, VL, VU, IL, IU, NS, SSAV, U, LDU, VT, LDVT, WORK, LWORK, RWORK, IWORK, IINFO )

                     // Compare U

                     DIF = ZERO
                     if ( M.GT.0 .AND. N.GT.0 ) {
                        if ( IJU.EQ.1 ) {
                           CALL CUNT03( 'C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, RWORK, DIF, IINFO )
                        }
                     }
                     RESULT( 27 ) = MAX( RESULT( 27 ), DIF )

                     // Compare VT

                     DIF = ZERO
                     if ( M.GT.0 .AND. N.GT.0 ) {
                        if ( IJVT.EQ.1 ) {
                           CALL CUNT03( 'R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, RWORK, DIF, IINFO )
                        }
                     }
                     RESULT( 28 ) = MAX( RESULT( 28 ), DIF )

                     // Compare S

                     DIF = ZERO
                     DIV = MAX( REAL( MNMIN )*ULP*S( 1 ), SLAMCH( 'Safe minimum' ) )
                     DO 150 I = 1, MNMIN - 1
                        IF( SSAV( I ).LT.SSAV( I+1 ) ) DIF = ULPINV                         IF( SSAV( I ).LT.ZERO ) DIF = ULPINV
                        DIF = MAX( DIF, ABS( SSAV( I )-S( I ) ) / DIV )
  150                CONTINUE
                     RESULT( 29) = MAX( RESULT( 29 ), DIF )
  160             CONTINUE
  170          CONTINUE

               // Do tests 8--10

               DO 180 I = 1, 4
                  ISEED2( I ) = ISEED( I )
  180          CONTINUE
               if ( MNMIN.LE.1 ) {
                  IL = 1
                  IU = MAX( 1, MNMIN )
               } else {
                  IL = 1 + INT( ( MNMIN-1 )*SLARND( 1, ISEED2 ) )
                  IU = 1 + INT( ( MNMIN-1 )*SLARND( 1, ISEED2 ) )
                  if ( IU.LT.IL ) {
                     ITEMP = IU
                     IU = IL
                     IL = ITEMP
                  }
               }
               CALL CLACPY( 'F', M, N, ASAV, LDA, A, LDA )
               SRNAMT = 'CGESVDX'
               CALL CGESVDX( 'V', 'V', 'I', M, N, A, LDA, VL, VU, IL, IU, NSI, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, IWORK, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9995 )'GESVDX', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               }

               RESULT( 30 ) = ZERO
               RESULT( 31 ) = ZERO
               RESULT( 32 ) = ZERO
               CALL CBDT05( M, N, ASAV, LDA, S, NSI, U, LDU, VT, LDVT, WORK, RESULT( 30 ) )
               if ( M.NE.0 .AND. N.NE.0 ) {
                  CALL CUNT01( 'Columns', M, NSI, U, LDU, WORK, LWORK, RWORK, RESULT( 31 ) )                   CALL CUNT01( 'Rows', NSI, N, VT, LDVT, WORK, LWORK, RWORK, RESULT( 32 ) )
               }

               // Do tests 11--13

               if ( MNMIN.GT.0 .AND. NSI.GT.1 ) {
                  if ( IL.NE.1 ) {
                     VU = SSAV( IL ) + MAX( HALF*ABS( SSAV( IL )-SSAV( IL-1 ) ), ULP*ANORM, TWO*RTUNFL )
                  } else {
                     VU = SSAV( 1 ) + MAX( HALF*ABS( SSAV( NS )-SSAV( 1 ) ), ULP*ANORM, TWO*RTUNFL )
                  }
                  if ( IU.NE.NS ) {
                     VL = SSAV( IU ) - MAX( ULP*ANORM, TWO*RTUNFL, HALF*ABS( SSAV( IU+1 )-SSAV( IU ) ) )
                  } else {
                     VL = SSAV( NS ) - MAX( ULP*ANORM, TWO*RTUNFL, HALF*ABS( SSAV( NS )-SSAV( 1 ) ) )
                  }
                  VL = MAX( VL,ZERO )
                  VU = MAX( VU,ZERO )
                  IF( VL.GE.VU ) VU = MAX( VU*2, VU+VL+HALF )
               } else {
                  VL = ZERO
                  VU = ONE
               }
               CALL CLACPY( 'F', M, N, ASAV, LDA, A, LDA )
               SRNAMT = 'CGESVDX'
               CALL CGESVDX( 'V', 'V', 'V', M, N, A, LDA, VL, VU, IL, IU, NSV, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, IWORK, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUNIT, FMT = 9995 )'GESVDX', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               }

               RESULT( 33 ) = ZERO
               RESULT( 34 ) = ZERO
               RESULT( 35 ) = ZERO
               CALL CBDT05( M, N, ASAV, LDA, S, NSV, U, LDU, VT, LDVT, WORK, RESULT( 33 ) )
               if ( M.NE.0 .AND. N.NE.0 ) {
                  CALL CUNT01( 'Columns', M, NSV, U, LDU, WORK, LWORK, RWORK, RESULT( 34 ) )                   CALL CUNT01( 'Rows', NSV, N, VT, LDVT, WORK, LWORK, RWORK, RESULT( 35 ) )
               }

               // End of Loop -- Check for RESULT(j) > THRESH

               NTEST = 0
               NFAIL = 0
               DO 190 J = 1, 39
                  IF( RESULT( J ).GE.ZERO ) NTEST = NTEST + 1                   IF( RESULT( J ).GE.THRESH ) NFAIL = NFAIL + 1
  190          CONTINUE

               IF( NFAIL.GT.0 ) NTESTF = NTESTF + 1
               if ( NTESTF.EQ.1 ) {
                  WRITE( NOUNIT, FMT = 9999 )
                  WRITE( NOUNIT, FMT = 9998 )THRESH
                  NTESTF = 2
               }

               DO 200 J = 1, 39
                  if ( RESULT( J ).GE.THRESH ) {
                     WRITE( NOUNIT, FMT = 9997 )M, N, JTYPE, IWSPC, IOLDSD, J, RESULT( J )
                  }
  200          CONTINUE

               NERRS = NERRS + NFAIL
               NTESTT = NTESTT + NTEST

  290       CONTINUE

  300    CONTINUE
  310 CONTINUE

      // Summary

      CALL ALASVM( 'CBD', NOUNIT, NERRS, NTESTT, 0 )

 9999 FORMAT( ' SVD -- Complex Singular Value Decomposition Driver ', / ' Matrix types (see CDRVBD for details):', / / ' 1 = Zero matrix', / ' 2 = Identity matrix', / ' 3 = Evenly spaced singular values near 1', / ' 4 = Evenly spaced singular values near underflow', / ' 5 = Evenly spaced singular values near overflow', / / ' Tests performed: ( A is dense, U and V are unitary,', / 19X, ' S is an array, and Upartial, VTpartial, and', / 19X, ' Spartial are partially computed U, VT and S),', / )
 9998 FORMAT( ' Tests performed with Test Threshold = ', F8.2, / ' CGESVD: ', / ' 1 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ', / ' 2 = | I - U**T U | / ( M ulp ) ', / ' 3 = | I - VT VT**T | / ( N ulp ) ', / ' 4 = 0 if S contains min(M,N) nonnegative values in', ' decreasing order, else 1/ulp', / ' 5 = | U - Upartial | / ( M ulp )', / ' 6 = | VT - VTpartial | / ( N ulp )', / ' 7 = | S - Spartial | / ( min(M,N) ulp |S| )', / ' CGESDD: ', / ' 8 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ', / ' 9 = | I - U**T U | / ( M ulp ) ', / '10 = | I - VT VT**T | / ( N ulp ) ', / '11 = 0 if S contains min(M,N) nonnegative values in', ' decreasing order, else 1/ulp', / '12 = | U - Upartial | / ( M ulp )', / '13 = | VT - VTpartial | / ( N ulp )', / '14 = | S - Spartial | / ( min(M,N) ulp |S| )', / ' CGESVJ: ', / / '15 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ', / '16 = | I - U**T U | / ( M ulp ) ', / '17 = | I - VT VT**T | / ( N ulp ) ', / '18 = 0 if S contains min(M,N) nonnegative values in', ' decreasing order, else 1/ulp', / ' CGESJV: ', / / '19 = | A - U diag(S) VT | / ( |A| max(M,N) ulp )',/ '20 = | I - U**T U | / ( M ulp ) ',/ '21 = | I - VT VT**T | / ( N ulp ) ',/ '22 = 0 if S contains min(M,N) nonnegative values in',' decreasing order, else 1/ulp',/ ' CGESVDX(V,V,A): ', /  '23 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ',/ '24 = | I - U**T U | / ( M ulp ) ',/ '25 = | I - VT VT**T | / ( N ulp ) ',/ '26 = 0 if S contains min(M,N) nonnegative values in',' decreasing order, else 1/ulp',/ '27 = | U - Upartial | / ( M ulp )',/ '28 = | VT - VTpartial | / ( N ulp )',/ '29 = | S - Spartial | / ( min(M,N) ulp |S| )',/ ' CGESVDX(V,V,I): ',/ '30 = | U**T A VT**T - diag(S) | / ( |A| max(M,N) ulp )',/ '31 = | I - U**T U | / ( M ulp ) ',/ '32 = | I - VT VT**T | / ( N ulp ) ',/ ' CGESVDX(V,V,V) ',/ '33 = | U**T A VT**T - diag(S) | / ( |A| max(M,N) ulp )',/ '34 = | I - U**T U | / ( M ulp ) ',/ '35 = | I - VT VT**T | / ( N ulp ) ',' CGESVDQ(H,N,N,A,A',/ '36 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ',/ '37 = | I - U**T U | / ( M ulp ) ',/ '38 = | I - VT VT**T | / ( N ulp ) ',/ '39 = 0 if S contains min(M,N) nonnegative values in',' decreasing order, else 1/ulp',/ / )
 9997 FORMAT( ' M=', I5, ', N=', I5, ', type ', I1, ', IWS=', I1, ', seed=', 4( I4, ',' ), ' test(', I2, ')=', G11.4 )
 9996 FORMAT( ' CDRVBD: ', A, ' returned INFO=', I6, '.', / 9X, 'M=', I6, ', N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )
 9995 FORMAT( ' CDRVBD: ', A, ' returned INFO=', I6, '.', / 9X, 'M=', I6, ', N=', I6, ', JTYPE=', I6, ', LSWORK=', I6, / 9X, 'ISEED=(', 3( I5, ',' ), I5, ')' )

      RETURN

      // End of CDRVBD

      }
