      SUBROUTINE SDRVBD( NSIZES, MM, NN, NTYPES, DOTYPE, ISEED, THRESH, A, LDA, U, LDU, VT, LDVT, ASAV, USAV, VTSAV, S, SSAV, E, WORK, LWORK, IWORK, NOUT, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      IMPLICIT NONE

      // .. Scalar Arguments ..
      int                INFO, LDA, LDU, LDVT, LWORK, NOUT, NSIZES, NTYPES;
      REAL               THRESH
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                ISEED( 4 ), IWORK( * ), MM( * ), NN( * );
      REAL               A( LDA, * ), ASAV( LDA, * ), E( * ), S( * ), SSAV( * ), U( LDU, * ), USAV( LDU, * ), VT( LDVT, * ), VTSAV( LDVT, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL              ZERO, ONE, TWO, HALF
      const              ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0, HALF = 0.5E0 ;
      int                MAXTYP;
      const              MAXTYP = 5 ;
      // ..
      // .. Local Scalars ..
      bool               BADMM, BADNN;
      String             JOBQ, JOBU, JOBVT, RANGE;
      String             PATH;
      int                I, IINFO, IJQ, IJU, IJVT, IL,IU, IWS, IWTMP, ITEMP, J, JSIZE, JTYPE, LSWORK, M, MINWRK, MMAX, MNMAX, MNMIN, MTYPES, N, NFAIL, NMAX, NS, NSI, NSV, NTEST;
      REAL               ANORM, DIF, DIV, OVFL, RTUNFL, ULP, ULPINV, UNFL, VL, VU
      // ..
      // .. Local Scalars for DGESVDQ ..
      int                LIWORK, LRWORK, NUMRANK;
      // ..
      // .. Local Arrays for DGESVDQ ..
      REAL               RWORK( 2 )
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
      // EXTERNAL ALASVM, SBDT01, SGEJSV, SGESDD, SGESVD, SGESVDQ, SGESVDX, SGESVJ, SLACPY, SLASET, SLATMS, SORT01, SORT03, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, REAL, INT, MAX, MIN
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NUNIT;
      // ..
      // .. Common blocks ..
      COMMON             / INFOC / INFOT, NUNIT, OK, LERR
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
         MINWRK = MAX( MINWRK, MAX( 3*MIN( MM( J ), NN( J ) )+MAX( MM( J ), NN( J ) ), 5*MIN( MM( J ), NN( J )-4 ) )+2*MIN( MM( J ), NN( J ) )**2 )
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
         xerbla('SDRVBD', -INFO );
         RETURN
      }

      // Initialize constants

      PATH( 1: 1 ) = 'Single precision'
      PATH( 2: 3 ) = 'BD'
      NFAIL = 0
      NTEST = 0
      UNFL = SLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      ULP = SLAMCH( 'Precision' )
      RTUNFL = SQRT( UNFL )
      ULPINV = ONE / ULP
      INFOT = 0

      // Loop over sizes, types

      DO 240 JSIZE = 1, NSIZES
         M = MM( JSIZE )
         N = NN( JSIZE )
         MNMIN = MIN( M, N )

         if ( NSIZES.NE.1 ) {
            MTYPES = MIN( MAXTYP, NTYPES )
         } else {
            MTYPES = MIN( MAXTYP+1, NTYPES )
         }

         DO 230 JTYPE = 1, MTYPES
            IF( .NOT.DOTYPE( JTYPE ) ) GO TO 230

            DO 20 J = 1, 4
               IOLDSD( J ) = ISEED( J )
   20       CONTINUE

            // Compute "A"

            IF( MTYPES.GT.MAXTYP ) GO TO 30

            if ( JTYPE.EQ.1 ) {

               // Zero matrix

               slaset('Full', M, N, ZERO, ZERO, A, LDA );

            } else if ( JTYPE.EQ.2 ) {

               // Identity matrix

               slaset('Full', M, N, ZERO, ONE, A, LDA );

            } else {

               // (Scaled) random matrix

               IF( JTYPE.EQ.3 ) ANORM = ONE                IF( JTYPE.EQ.4 ) ANORM = UNFL / ULP                IF( JTYPE.EQ.5 ) ANORM = OVFL*ULP                CALL SLATMS( M, N, 'U', ISEED, 'N', S, 4, REAL( MNMIN ), ANORM, M-1, N-1, 'N', A, LDA, WORK, IINFO )
               if ( IINFO.NE.0 ) {
                  WRITE( NOUT, FMT = 9996 )'Generator', IINFO, M, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               }
            }

   30       CONTINUE
            slacpy('F', M, N, A, LDA, ASAV, LDA );

            // Do for minimal and adequate (for blocking) workspace

            DO 220 IWS = 1, 4

               DO 40 J = 1, 32
                  RESULT( J ) = -ONE
   40          CONTINUE

               // Test SGESVD: Factorize A

               IWTMP = MAX( 3*MIN( M, N )+MAX( M, N ), 5*MIN( M, N ) )
               LSWORK = IWTMP + ( IWS-1 )*( LWORK-IWTMP ) / 3
               LSWORK = MIN( LSWORK, LWORK )
               LSWORK = MAX( LSWORK, 1 )
               IF( IWS.EQ.4 ) LSWORK = LWORK

               IF( IWS.GT.1 ) CALL SLACPY( 'F', M, N, ASAV, LDA, A, LDA )
               SRNAMT = 'SGESVD'
               sgesvd('A', 'A', M, N, A, LDA, SSAV, USAV, LDU, VTSAV, LDVT, WORK, LSWORK, IINFO );
               if ( IINFO.NE.0 ) {
                  WRITE( NOUT, FMT = 9995 )'GESVD', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               }

               // Do tests 1--4

               sbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RESULT( 1 ) );
               if ( M.NE.0 .AND. N.NE.0 ) {
                  sort01('Columns', M, M, USAV, LDU, WORK, LWORK, RESULT( 2 ) )                   CALL SORT01( 'Rows', N, N, VTSAV, LDVT, WORK, LWORK, RESULT( 3 ) );
               }
               RESULT( 4 ) = ZERO
               DO 50 I = 1, MNMIN - 1
                  IF( SSAV( I ).LT.SSAV( I+1 ) ) RESULT( 4 ) = ULPINV                   IF( SSAV( I ).LT.ZERO ) RESULT( 4 ) = ULPINV
   50          CONTINUE
               if ( MNMIN.GE.1 ) {
                  IF( SSAV( MNMIN ).LT.ZERO ) RESULT( 4 ) = ULPINV
               }

               // Do partial SVDs, comparing to SSAV, USAV, and VTSAV

               RESULT( 5 ) = ZERO
               RESULT( 6 ) = ZERO
               RESULT( 7 ) = ZERO
               DO 80 IJU = 0, 3
                  DO 70 IJVT = 0, 3
                     IF( ( IJU.EQ.3 .AND. IJVT.EQ.3 ) .OR. ( IJU.EQ.1 .AND. IJVT.EQ.1 ) )GO TO 70
                     JOBU = CJOB( IJU+1 )
                     JOBVT = CJOB( IJVT+1 )
                     slacpy('F', M, N, ASAV, LDA, A, LDA );
                     SRNAMT = 'SGESVD'
                     sgesvd(JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LSWORK, IINFO );

                     // Compare U

                     DIF = ZERO
                     if ( M.GT.0 .AND. N.GT.0 ) {
                        if ( IJU.EQ.1 ) {
                           sort03('C', M, MNMIN, M, MNMIN, USAV, LDU, A, LDA, WORK, LWORK, DIF, IINFO );
                        } else if ( IJU.EQ.2 ) {
                           sort03('C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, DIF, IINFO );
                        } else if ( IJU.EQ.3 ) {
                           sort03('C', M, M, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, DIF, IINFO );
                        }
                     }
                     RESULT( 5 ) = MAX( RESULT( 5 ), DIF )

                     // Compare VT

                     DIF = ZERO
                     if ( M.GT.0 .AND. N.GT.0 ) {
                        if ( IJVT.EQ.1 ) {
                           sort03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, A, LDA, WORK, LWORK, DIF, IINFO );
                        } else if ( IJVT.EQ.2 ) {
                           sort03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, DIF, IINFO );
                        } else if ( IJVT.EQ.3 ) {
                           sort03('R', N, N, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, DIF, IINFO );
                        }
                     }
                     RESULT( 6 ) = MAX( RESULT( 6 ), DIF )

                     // Compare S

                     DIF = ZERO
                     DIV = MAX( MNMIN*ULP*S( 1 ), UNFL )
                     DO 60 I = 1, MNMIN - 1
                        IF( SSAV( I ).LT.SSAV( I+1 ) ) DIF = ULPINV                         IF( SSAV( I ).LT.ZERO ) DIF = ULPINV
                        DIF = MAX( DIF, ABS( SSAV( I )-S( I ) ) / DIV )
   60                CONTINUE
                     RESULT( 7 ) = MAX( RESULT( 7 ), DIF )
   70             CONTINUE
   80          CONTINUE

               // Test SGESDD: Factorize A

               IWTMP = 5*MNMIN*MNMIN + 9*MNMIN + MAX( M, N )
               LSWORK = IWTMP + ( IWS-1 )*( LWORK-IWTMP ) / 3
               LSWORK = MIN( LSWORK, LWORK )
               LSWORK = MAX( LSWORK, 1 )
               IF( IWS.EQ.4 ) LSWORK = LWORK

               slacpy('F', M, N, ASAV, LDA, A, LDA );
               SRNAMT = 'SGESDD'
               sgesdd('A', M, N, A, LDA, SSAV, USAV, LDU, VTSAV, LDVT, WORK, LSWORK, IWORK, IINFO );
               if ( IINFO.NE.0 ) {
                  WRITE( NOUT, FMT = 9995 )'GESDD', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               }

               // Do tests 8--11

               sbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RESULT( 8 ) );
               if ( M.NE.0 .AND. N.NE.0 ) {
                  sort01('Columns', M, M, USAV, LDU, WORK, LWORK, RESULT( 9 ) )                   CALL SORT01( 'Rows', N, N, VTSAV, LDVT, WORK, LWORK, RESULT( 10 ) );
               }
               RESULT( 11 ) = ZERO
               DO 90 I = 1, MNMIN - 1
                  IF( SSAV( I ).LT.SSAV( I+1 ) ) RESULT( 11 ) = ULPINV                   IF( SSAV( I ).LT.ZERO ) RESULT( 11 ) = ULPINV
   90          CONTINUE
               if ( MNMIN.GE.1 ) {
                  IF( SSAV( MNMIN ).LT.ZERO ) RESULT( 11 ) = ULPINV
               }

               // Do partial SVDs, comparing to SSAV, USAV, and VTSAV

               RESULT( 12 ) = ZERO
               RESULT( 13 ) = ZERO
               RESULT( 14 ) = ZERO
               DO 110 IJQ = 0, 2
                  JOBQ = CJOB( IJQ+1 )
                  slacpy('F', M, N, ASAV, LDA, A, LDA );
                  SRNAMT = 'SGESDD'
                  sgesdd(JOBQ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LSWORK, IWORK, IINFO );

                  // Compare U

                  DIF = ZERO
                  if ( M.GT.0 .AND. N.GT.0 ) {
                     if ( IJQ.EQ.1 ) {
                        if ( M.GE.N ) {
                           sort03('C', M, MNMIN, M, MNMIN, USAV, LDU, A, LDA, WORK, LWORK, DIF, INFO );
                        } else {
                           sort03('C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, DIF, INFO );
                        }
                     } else if ( IJQ.EQ.2 ) {
                        sort03('C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, DIF, INFO );
                     }
                  }
                  RESULT( 12 ) = MAX( RESULT( 12 ), DIF )

                  // Compare VT

                  DIF = ZERO
                  if ( M.GT.0 .AND. N.GT.0 ) {
                     if ( IJQ.EQ.1 ) {
                        if ( M.GE.N ) {
                           sort03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, DIF, INFO );
                        } else {
                           sort03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, A, LDA, WORK, LWORK, DIF, INFO );
                        }
                     } else if ( IJQ.EQ.2 ) {
                        sort03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, DIF, INFO );
                     }
                  }
                  RESULT( 13 ) = MAX( RESULT( 13 ), DIF )

                  // Compare S

                  DIF = ZERO
                  DIV = MAX( MNMIN*ULP*S( 1 ), UNFL )
                  DO 100 I = 1, MNMIN - 1
                     IF( SSAV( I ).LT.SSAV( I+1 ) ) DIF = ULPINV                      IF( SSAV( I ).LT.ZERO ) DIF = ULPINV
                     DIF = MAX( DIF, ABS( SSAV( I )-S( I ) ) / DIV )
  100             CONTINUE
                  RESULT( 14 ) = MAX( RESULT( 14 ), DIF )
  110          CONTINUE

               // Test SGESVDQ
               // Note: SGESVDQ only works for M >= N

               RESULT( 36 ) = ZERO
               RESULT( 37 ) = ZERO
               RESULT( 38 ) = ZERO
               RESULT( 39 ) = ZERO

               if ( M.GE.N ) {
                  IWTMP = 5*MNMIN*MNMIN + 9*MNMIN + MAX( M, N )
                  LSWORK = IWTMP + ( IWS-1 )*( LWORK-IWTMP ) / 3
                  LSWORK = MIN( LSWORK, LWORK )
                  LSWORK = MAX( LSWORK, 1 )
                  IF( IWS.EQ.4 ) LSWORK = LWORK

                  slacpy('F', M, N, ASAV, LDA, A, LDA );
                  SRNAMT = 'SGESVDQ'

                  LRWORK = 2
                  LIWORK = MAX( N, 1 )
                  sgesvdq('H', 'N', 'N', 'A', 'A',  M, N, A, LDA, SSAV, USAV, LDU, VTSAV, LDVT, NUMRANK, IWORK, LIWORK, WORK, LWORK, RWORK, LRWORK, IINFO );

                  if ( IINFO.NE.0 ) {
                     WRITE( NOUT, FMT = 9995 )'SGESVDQ', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                     INFO = ABS( IINFO )
                     RETURN
                  }

                  // Do tests 36--39

                  sbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RESULT( 36 ) );
                  if ( M.NE.0 .AND. N.NE.0 ) {
                     sort01('Columns', M, M, USAV, LDU, WORK, LWORK, RESULT( 37 ) )                      CALL SORT01( 'Rows', N, N, VTSAV, LDVT, WORK, LWORK, RESULT( 38 ) );
                  }
                  RESULT( 39 ) = ZERO
                  DO 199 I = 1, MNMIN - 1
                     IF( SSAV( I ).LT.SSAV( I+1 ) ) RESULT( 39 ) = ULPINV                      IF( SSAV( I ).LT.ZERO ) RESULT( 39 ) = ULPINV
  199             CONTINUE
                  if ( MNMIN.GE.1 ) {
                     IF( SSAV( MNMIN ).LT.ZERO ) RESULT( 39 ) = ULPINV
                  }
               }

               // Test SGESVJ
               // Note: SGESVJ only works for M >= N

               RESULT( 15 ) = ZERO
               RESULT( 16 ) = ZERO
               RESULT( 17 ) = ZERO
               RESULT( 18 ) = ZERO

               if ( M.GE.N ) {
                  IWTMP = 5*MNMIN*MNMIN + 9*MNMIN + MAX( M, N )
                  LSWORK = IWTMP + ( IWS-1 )*( LWORK-IWTMP ) / 3
                  LSWORK = MIN( LSWORK, LWORK )
                  LSWORK = MAX( LSWORK, 1 )
                  IF( IWS.EQ.4 ) LSWORK = LWORK

                  slacpy('F', M, N, ASAV, LDA, USAV, LDA );
                  SRNAMT = 'SGESVJ'
                  sgesvj('G', 'U', 'V', M, N, USAV, LDA, SSAV, 0, A, LDVT, WORK, LWORK, INFO );

                  // SGESVJ returns V not VT

                  DO J=1,N
                     DO I=1,N
                        VTSAV(J,I) = A(I,J)
                     END DO
                  END DO

                  if ( IINFO.NE.0 ) {
                     WRITE( NOUT, FMT = 9995 )'GESVJ', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                     INFO = ABS( IINFO )
                     RETURN
                  }

                  // Do tests 15--18

                  sbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RESULT( 15 ) );
                  if ( M.NE.0 .AND. N.NE.0 ) {
                     sort01('Columns', M, M, USAV, LDU, WORK, LWORK, RESULT( 16 ) )                      CALL SORT01( 'Rows', N, N, VTSAV, LDVT, WORK, LWORK, RESULT( 17 ) );
                  }
                  RESULT( 18 ) = ZERO
                  DO 120 I = 1, MNMIN - 1
                     IF( SSAV( I ).LT.SSAV( I+1 ) ) RESULT( 18 ) = ULPINV                      IF( SSAV( I ).LT.ZERO ) RESULT( 18 ) = ULPINV
  120             CONTINUE
                  if ( MNMIN.GE.1 ) {
                     IF( SSAV( MNMIN ).LT.ZERO ) RESULT( 18 ) = ULPINV
                  }
               }

               // Test SGEJSV
               // Note: SGEJSV only works for M >= N

               RESULT( 19 ) = ZERO
               RESULT( 20 ) = ZERO
               RESULT( 21 ) = ZERO
               RESULT( 22 ) = ZERO
               if ( M.GE.N ) {
                  IWTMP = 5*MNMIN*MNMIN + 9*MNMIN + MAX( M, N )
                  LSWORK = IWTMP + ( IWS-1 )*( LWORK-IWTMP ) / 3
                  LSWORK = MIN( LSWORK, LWORK )
                  LSWORK = MAX( LSWORK, 1 )
                  IF( IWS.EQ.4 ) LSWORK = LWORK

                  slacpy('F', M, N, ASAV, LDA, VTSAV, LDA );
                  SRNAMT = 'SGEJSV'
                  sgejsv('G', 'U', 'V', 'R', 'N', 'N', M, N, VTSAV, LDA, SSAV, USAV, LDU, A, LDVT, WORK, LWORK, IWORK, INFO );

                  // SGEJSV returns V not VT

                  DO 140 J=1,N
                     DO 130 I=1,N
                        VTSAV(J,I) = A(I,J)
  130                END DO
  140             END DO

                  if ( IINFO.NE.0 ) {
                     WRITE( NOUT, FMT = 9995 )'GEJSV', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                     INFO = ABS( IINFO )
                     RETURN
                  }

                  // Do tests 19--22

                  sbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RESULT( 19 ) );
                  if ( M.NE.0 .AND. N.NE.0 ) {
                     sort01('Columns', M, M, USAV, LDU, WORK, LWORK, RESULT( 20 ) )                      CALL SORT01( 'Rows', N, N, VTSAV, LDVT, WORK, LWORK, RESULT( 21 ) );
                  }
                  RESULT( 22 ) = ZERO
                  DO 150 I = 1, MNMIN - 1
                     IF( SSAV( I ).LT.SSAV( I+1 ) ) RESULT( 22 ) = ULPINV                      IF( SSAV( I ).LT.ZERO ) RESULT( 22 ) = ULPINV
  150             CONTINUE
                  if ( MNMIN.GE.1 ) {
                     IF( SSAV( MNMIN ).LT.ZERO ) RESULT( 22 ) = ULPINV
                  }
               }

               // Test SGESVDX

               slacpy('F', M, N, ASAV, LDA, A, LDA );
               sgesvdx('V', 'V', 'A', M, N, A, LDA, VL, VU, IL, IU, NS, SSAV, USAV, LDU, VTSAV, LDVT, WORK, LWORK, IWORK, IINFO );
               if ( IINFO.NE.0 ) {
                  WRITE( NOUT, FMT = 9995 )'GESVDX', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               }

               // Do tests 23--29

               RESULT( 23 ) = ZERO
               RESULT( 24 ) = ZERO
               RESULT( 25 ) = ZERO
               sbdt01(M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RESULT( 23 ) );
               if ( M.NE.0 .AND. N.NE.0 ) {
                  sort01('Columns', M, M, USAV, LDU, WORK, LWORK, RESULT( 24 ) )                   CALL SORT01( 'Rows', N, N, VTSAV, LDVT, WORK, LWORK, RESULT( 25 ) );
               }
               RESULT( 26 ) = ZERO
               DO 160 I = 1, MNMIN - 1
                  IF( SSAV( I ).LT.SSAV( I+1 ) ) RESULT( 26 ) = ULPINV                   IF( SSAV( I ).LT.ZERO ) RESULT( 26 ) = ULPINV
  160          CONTINUE
               if ( MNMIN.GE.1 ) {
                  IF( SSAV( MNMIN ).LT.ZERO ) RESULT( 26 ) = ULPINV
               }

               // Do partial SVDs, comparing to SSAV, USAV, and VTSAV

               RESULT( 27 ) = ZERO
               RESULT( 28 ) = ZERO
               RESULT( 29 ) = ZERO
               DO 180 IJU = 0, 1
                  DO 170 IJVT = 0, 1
                     IF( ( IJU.EQ.0 .AND. IJVT.EQ.0 ) .OR. ( IJU.EQ.1 .AND. IJVT.EQ.1 ) )GO TO 170
                     JOBU = CJOBV( IJU+1 )
                     JOBVT = CJOBV( IJVT+1 )
                     RANGE = CJOBR( 1 )
                     slacpy('F', M, N, ASAV, LDA, A, LDA );
                     sgesvdx(JOBU, JOBVT, RANGE, M, N, A, LDA, VL, VU, IL, IU, NS, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, IINFO );

                     // Compare U

                     DIF = ZERO
                     if ( M.GT.0 .AND. N.GT.0 ) {
                        if ( IJU.EQ.1 ) {
                           sort03('C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, DIF, IINFO );
                        }
                     }
                     RESULT( 27 ) = MAX( RESULT( 27 ), DIF )

                     // Compare VT

                     DIF = ZERO
                     if ( M.GT.0 .AND. N.GT.0 ) {
                        if ( IJVT.EQ.1 ) {
                           sort03('R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, DIF, IINFO );
                        }
                     }
                     RESULT( 28 ) = MAX( RESULT( 28 ), DIF )

                     // Compare S

                     DIF = ZERO
                     DIV = MAX( MNMIN*ULP*S( 1 ), UNFL )
                     DO 190 I = 1, MNMIN - 1
                        IF( SSAV( I ).LT.SSAV( I+1 ) ) DIF = ULPINV                         IF( SSAV( I ).LT.ZERO ) DIF = ULPINV
                        DIF = MAX( DIF, ABS( SSAV( I )-S( I ) ) / DIV )
  190                CONTINUE
                     RESULT( 29 ) = MAX( RESULT( 29 ), DIF )
  170             CONTINUE
  180          CONTINUE

               // Do tests 30--32: SGESVDX( 'V', 'V', 'I' )

               DO 200 I = 1, 4
                  ISEED2( I ) = ISEED( I )
  200          CONTINUE
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
               slacpy('F', M, N, ASAV, LDA, A, LDA );
               sgesvdx('V', 'V', 'I', M, N, A, LDA, VL, VU, IL, IU, NSI, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, IINFO );
               if ( IINFO.NE.0 ) {
                  WRITE( NOUT, FMT = 9995 )'GESVDX', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               }

               RESULT( 30 ) = ZERO
               RESULT( 31 ) = ZERO
               RESULT( 32 ) = ZERO
               sbdt05(M, N, ASAV, LDA, S, NSI, U, LDU, VT, LDVT, WORK, RESULT( 30 ) )                CALL SORT01( 'Columns', M, NSI, U, LDU, WORK, LWORK, RESULT( 31 ) )                CALL SORT01( 'Rows', NSI, N, VT, LDVT, WORK, LWORK, RESULT( 32 ) );

               // Do tests 33--35: SGESVDX( 'V', 'V', 'V' )

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
               slacpy('F', M, N, ASAV, LDA, A, LDA );
               sgesvdx('V', 'V', 'V', M, N, A, LDA, VL, VU, IL, IU, NSV, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, IINFO );
               if ( IINFO.NE.0 ) {
                  WRITE( NOUT, FMT = 9995 )'GESVDX', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               }

               RESULT( 33 ) = ZERO
               RESULT( 34 ) = ZERO
               RESULT( 35 ) = ZERO
               sbdt05(M, N, ASAV, LDA, S, NSV, U, LDU, VT, LDVT, WORK, RESULT( 33 ) )                CALL SORT01( 'Columns', M, NSV, U, LDU, WORK, LWORK, RESULT( 34 ) )                CALL SORT01( 'Rows', NSV, N, VT, LDVT, WORK, LWORK, RESULT( 35 ) );

               // End of Loop -- Check for RESULT(j) > THRESH

               DO 210 J = 1, 39
                  if ( RESULT( J ).GE.THRESH ) {
                     if ( NFAIL.EQ.0 ) {
                        WRITE( NOUT, FMT = 9999 )
                        WRITE( NOUT, FMT = 9998 )
                     }
                     WRITE( NOUT, FMT = 9997 )M, N, JTYPE, IWS, IOLDSD, J, RESULT( J )
                     NFAIL = NFAIL + 1
                  }
  210          CONTINUE
               NTEST = NTEST + 39
  220       CONTINUE
  230    CONTINUE
  240 CONTINUE

      // Summary

      alasvm(PATH, NOUT, NFAIL, NTEST, 0 );

 9999 FORMAT( ' SVD -- Real Singular Value Decomposition Driver ', / ' Matrix types (see SDRVBD for details):', / / ' 1 = Zero matrix', / ' 2 = Identity matrix', / ' 3 = Evenly spaced singular values near 1', / ' 4 = Evenly spaced singular values near underflow', / ' 5 = Evenly spaced singular values near overflow', / / ' Tests performed: ( A is dense, U and V are orthogonal,', / 19X, ' S is an array, and Upartial, VTpartial, and', / 19X, ' Spartial are partially computed U, VT and S),', / )
 9998 FORMAT( ' 1 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ', / ' 2 = | I - U**T U | / ( M ulp ) ', / ' 3 = | I - VT VT**T | / ( N ulp ) ', / ' 4 = 0 if S contains min(M,N) nonnegative values in', ' decreasing order, else 1/ulp', / ' 5 = | U - Upartial | / ( M ulp )', / ' 6 = | VT - VTpartial | / ( N ulp )', / ' 7 = | S - Spartial | / ( min(M,N) ulp |S| )', / ' 8 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ', / ' 9 = | I - U**T U | / ( M ulp ) ', / '10 = | I - VT VT**T | / ( N ulp ) ', / '11 = 0 if S contains min(M,N) nonnegative values in', ' decreasing order, else 1/ulp', / '12 = | U - Upartial | / ( M ulp )', / '13 = | VT - VTpartial | / ( N ulp )', / '14 = | S - Spartial | / ( min(M,N) ulp |S| )', / '15 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ', / '16 = | I - U**T U | / ( M ulp ) ', / '17 = | I - VT VT**T | / ( N ulp ) ', / '18 = 0 if S contains min(M,N) nonnegative values in', ' decreasing order, else 1/ulp', / '19 = | U - Upartial | / ( M ulp )', / '20 = | VT - VTpartial | / ( N ulp )',/ '21 = | S - Spartial | / ( min(M,N) ulp |S| )',/ '22 = 0 if S contains min(M,N) nonnegative values in',' decreasing order, else 1/ulp',' SGESVDX(V,V,A) ',/ '23 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ),'/ '24 = | I - U**T U | / ( M ulp ) ',/ '25 = | I - VT VT**T | / ( N ulp ) ',/ '26 = 0 if S contains min(M,N) nonnegative values in',' decreasing order, else 1/ulp',/ '27 = | U - Upartial | / ( M ulp )',/ '28 = | VT - VTpartial | / ( N ulp )',/ '29 = | S - Spartial | / ( min(M,N) ulp |S| )',/ '30 = | U**T A VT**T - diag(S) | / ( |A| max(M,N) ulp ),',' SGESVDX(V,V,I) ',/ '31 = | I - U**T U | / ( M ulp ) ',/ '32 = | I - VT VT**T | / ( N ulp ) ',/ '33 = | U**T A VT**T - diag(S) | / ( |A| max(M,N) ulp ),',' SGESVDX(V,V,V) ',/ '34 = | I - U**T U | / ( M ulp ) ',/ '35 = | I - VT VT**T | / ( N ulp ) ',' SGESVDQ(H,N,N,A,A',/ '36 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ',/ '37 = | I - U**T U | / ( M ulp ) ',/ '38 = | I - VT VT**T | / ( N ulp ) ',/ '39 = 0 if S contains min(M,N) nonnegative values in',' decreasing order, else 1/ulp',/ / )
 9997 FORMAT( ' M=', I5, ', N=', I5, ', type ', I1, ', IWS=', I1, ', seed=', 4( I4, ',' ), ' test(', I2, ')=', G11.4 )
 9996 FORMAT( ' SDRVBD: ', A, ' returned INFO=', I6, '.', / 9X, 'M=', I6, ', N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )
 9995 FORMAT( ' SDRVBD: ', A, ' returned INFO=', I6, '.', / 9X, 'M=', I6, ', N=', I6, ', JTYPE=', I6, ', LSWORK=', I6, / 9X, 'ISEED=(', 3( I5, ',' ), I5, ')' )

      RETURN

      // End of SDRVBD

      }
