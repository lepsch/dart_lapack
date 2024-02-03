      SUBROUTINE DDRVBD( NSIZES, MM, NN, NTYPES, DOTYPE, ISEED, THRESH, A, LDA, U, LDU, VT, LDVT, ASAV, USAV, VTSAV, S, SSAV, E, WORK, LWORK, IWORK, NOUT, INFO )
*
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                INFO, LDA, LDU, LDVT, LWORK, NOUT, NSIZES, NTYPES;
      double             THRESH;
*     ..
*     .. Array Arguments ..
      bool               DOTYPE( * );
      int                ISEED( 4 ), IWORK( * ), MM( * ), NN( * );
      double             A( LDA, * ), ASAV( LDA, * ), E( * ), S( * ), SSAV( * ), U( LDU, * ), USAV( LDU, * ), VT( LDVT, * ), VTSAV( LDVT, * ), WORK( * );
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      double            ZERO, ONE, TWO, HALF;
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, HALF = 0.5D0 )
      int                MAXTYP;
      PARAMETER          ( MAXTYP = 5 )
*     ..
*     .. Local Scalars ..
      bool               BADMM, BADNN;
      String             JOBQ, JOBU, JOBVT, RANGE;
      String             PATH;
      int                I, IINFO, IJQ, IJU, IJVT, IL,IU, IWS, IWTMP, ITEMP, J, JSIZE, JTYPE, LSWORK, M, MINWRK, MMAX, MNMAX, MNMIN, MTYPES, N, NFAIL, NMAX, NS, NSI, NSV, NTEST;
      double             ANORM, DIF, DIV, OVFL, RTUNFL, ULP, ULPINV, UNFL, VL, VU;
*     ..
*     .. Local Scalars for DGESVDQ ..
      int                LIWORK, LRWORK, NUMRANK;
*     ..
*     .. Local Arrays for DGESVDQ ..
      double             RWORK( 2 );
*     ..
*     .. Local Arrays ..
      String             CJOB( 4 ), CJOBR( 3 ), CJOBV( 2 );
      int                IOLDSD( 4 ), ISEED2( 4 );
      double             RESULT( 39 );
*     ..
*     .. External Functions ..
      double             DLAMCH, DLARND;
      EXTERNAL           DLAMCH, DLARND
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALASVM, DBDT01, DGEJSV, DGESDD, DGESVD, DGESVDQ, DGESVDX, DGESVJ, DLACPY, DLASET, DLATMS, DORT01, DORT03, XERBLA
*     ..
*     .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, INT, MAX, MIN
*     ..
*     .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NUNIT;
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, NUNIT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Data statements ..
      DATA               CJOB / 'N', 'O', 'S', 'A' /
      DATA               CJOBR / 'A', 'V', 'I' /
      DATA               CJOBV / 'N', 'V' /
*     ..
*     .. Executable Statements ..
*
*     Check for errors
*
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
*
*     Check for errors
*
      IF( NSIZES.LT.0 ) THEN
         INFO = -1
      ELSE IF( BADMM ) THEN
         INFO = -2
      ELSE IF( BADNN ) THEN
         INFO = -3
      ELSE IF( NTYPES.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, MMAX ) ) THEN
         INFO = -10
      ELSE IF( LDU.LT.MAX( 1, MMAX ) ) THEN
         INFO = -12
      ELSE IF( LDVT.LT.MAX( 1, NMAX ) ) THEN
         INFO = -14
      ELSE IF( MINWRK.GT.LWORK ) THEN
         INFO = -21
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DDRVBD', -INFO )
         RETURN
      END IF
*
*     Initialize constants
*
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
*
*     Loop over sizes, types
*
      DO 240 JSIZE = 1, NSIZES
         M = MM( JSIZE )
         N = NN( JSIZE )
         MNMIN = MIN( M, N )
*
         IF( NSIZES.NE.1 ) THEN
            MTYPES = MIN( MAXTYP, NTYPES )
         ELSE
            MTYPES = MIN( MAXTYP+1, NTYPES )
         END IF
*
         DO 230 JTYPE = 1, MTYPES
            IF( .NOT.DOTYPE( JTYPE ) ) GO TO 230
*
            DO 20 J = 1, 4
               IOLDSD( J ) = ISEED( J )
   20       CONTINUE
*
*           Compute "A"
*
            IF( MTYPES.GT.MAXTYP ) GO TO 30
*
            IF( JTYPE.EQ.1 ) THEN
*
*              Zero matrix
*
               CALL DLASET( 'Full', M, N, ZERO, ZERO, A, LDA )
*
            ELSE IF( JTYPE.EQ.2 ) THEN
*
*              Identity matrix
*
               CALL DLASET( 'Full', M, N, ZERO, ONE, A, LDA )
*
            ELSE
*
*              (Scaled) random matrix
*
               IF( JTYPE.EQ.3 ) ANORM = ONE                IF( JTYPE.EQ.4 ) ANORM = UNFL / ULP                IF( JTYPE.EQ.5 ) ANORM = OVFL*ULP                CALL DLATMS( M, N, 'U', ISEED, 'N', S, 4, DBLE( MNMIN ), ANORM, M-1, N-1, 'N', A, LDA, WORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUT, FMT = 9996 )'Generator', IINFO, M, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               END IF
            END IF
*
   30       CONTINUE
            CALL DLACPY( 'F', M, N, A, LDA, ASAV, LDA )
*
*           Do for minimal and adequate (for blocking) workspace
*
            DO 220 IWS = 1, 4
*
               DO 40 J = 1, 32
                  RESULT( J ) = -ONE
   40          CONTINUE
*
*              Test DGESVD: Factorize A
*
               IWTMP = MAX( 3*MIN( M, N )+MAX( M, N ), 5*MIN( M, N ) )
               LSWORK = IWTMP + ( IWS-1 )*( LWORK-IWTMP ) / 3
               LSWORK = MIN( LSWORK, LWORK )
               LSWORK = MAX( LSWORK, 1 )
               IF( IWS.EQ.4 ) LSWORK = LWORK
*
               IF( IWS.GT.1 ) CALL DLACPY( 'F', M, N, ASAV, LDA, A, LDA )
               SRNAMT = 'DGESVD'
               CALL DGESVD( 'A', 'A', M, N, A, LDA, SSAV, USAV, LDU, VTSAV, LDVT, WORK, LSWORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUT, FMT = 9995 )'GESVD', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               END IF
*
*              Do tests 1--4
*
               CALL DBDT01( M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RESULT( 1 ) )
               IF( M.NE.0 .AND. N.NE.0 ) THEN
                  CALL DORT01( 'Columns', M, M, USAV, LDU, WORK, LWORK, RESULT( 2 ) )                   CALL DORT01( 'Rows', N, N, VTSAV, LDVT, WORK, LWORK, RESULT( 3 ) )
               END IF
               RESULT( 4 ) = ZERO
               DO 50 I = 1, MNMIN - 1
                  IF( SSAV( I ).LT.SSAV( I+1 ) ) RESULT( 4 ) = ULPINV                   IF( SSAV( I ).LT.ZERO ) RESULT( 4 ) = ULPINV
   50          CONTINUE
               IF( MNMIN.GE.1 ) THEN
                  IF( SSAV( MNMIN ).LT.ZERO ) RESULT( 4 ) = ULPINV
               END IF
*
*              Do partial SVDs, comparing to SSAV, USAV, and VTSAV
*
               RESULT( 5 ) = ZERO
               RESULT( 6 ) = ZERO
               RESULT( 7 ) = ZERO
               DO 80 IJU = 0, 3
                  DO 70 IJVT = 0, 3
                     IF( ( IJU.EQ.3 .AND. IJVT.EQ.3 ) .OR. ( IJU.EQ.1 .AND. IJVT.EQ.1 ) )GO TO 70
                     JOBU = CJOB( IJU+1 )
                     JOBVT = CJOB( IJVT+1 )
                     CALL DLACPY( 'F', M, N, ASAV, LDA, A, LDA )
                     SRNAMT = 'DGESVD'
                     CALL DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LSWORK, IINFO )
*
*                    Compare U
*
                     DIF = ZERO
                     IF( M.GT.0 .AND. N.GT.0 ) THEN
                        IF( IJU.EQ.1 ) THEN
                           CALL DORT03( 'C', M, MNMIN, M, MNMIN, USAV, LDU, A, LDA, WORK, LWORK, DIF, IINFO )
                        ELSE IF( IJU.EQ.2 ) THEN
                           CALL DORT03( 'C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, DIF, IINFO )
                        ELSE IF( IJU.EQ.3 ) THEN
                           CALL DORT03( 'C', M, M, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, DIF, IINFO )
                        END IF
                     END IF
                     RESULT( 5 ) = MAX( RESULT( 5 ), DIF )
*
*                    Compare VT
*
                     DIF = ZERO
                     IF( M.GT.0 .AND. N.GT.0 ) THEN
                        IF( IJVT.EQ.1 ) THEN
                           CALL DORT03( 'R', N, MNMIN, N, MNMIN, VTSAV, LDVT, A, LDA, WORK, LWORK, DIF, IINFO )
                        ELSE IF( IJVT.EQ.2 ) THEN
                           CALL DORT03( 'R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, DIF, IINFO )
                        ELSE IF( IJVT.EQ.3 ) THEN
                           CALL DORT03( 'R', N, N, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, DIF, IINFO )
                        END IF
                     END IF
                     RESULT( 6 ) = MAX( RESULT( 6 ), DIF )
*
*                    Compare S
*
                     DIF = ZERO
                     DIV = MAX( MNMIN*ULP*S( 1 ), UNFL )
                     DO 60 I = 1, MNMIN - 1
                        IF( SSAV( I ).LT.SSAV( I+1 ) ) DIF = ULPINV                         IF( SSAV( I ).LT.ZERO ) DIF = ULPINV
                        DIF = MAX( DIF, ABS( SSAV( I )-S( I ) ) / DIV )
   60                CONTINUE
                     RESULT( 7 ) = MAX( RESULT( 7 ), DIF )
   70             CONTINUE
   80          CONTINUE
*
*              Test DGESDD: Factorize A
*
               IWTMP = 5*MNMIN*MNMIN + 9*MNMIN + MAX( M, N )
               LSWORK = IWTMP + ( IWS-1 )*( LWORK-IWTMP ) / 3
               LSWORK = MIN( LSWORK, LWORK )
               LSWORK = MAX( LSWORK, 1 )
               IF( IWS.EQ.4 ) LSWORK = LWORK
*
               CALL DLACPY( 'F', M, N, ASAV, LDA, A, LDA )
               SRNAMT = 'DGESDD'
               CALL DGESDD( 'A', M, N, A, LDA, SSAV, USAV, LDU, VTSAV, LDVT, WORK, LSWORK, IWORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUT, FMT = 9995 )'GESDD', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               END IF
*
*              Do tests 8--11
*
               CALL DBDT01( M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RESULT( 8 ) )
               IF( M.NE.0 .AND. N.NE.0 ) THEN
                  CALL DORT01( 'Columns', M, M, USAV, LDU, WORK, LWORK, RESULT( 9 ) )                   CALL DORT01( 'Rows', N, N, VTSAV, LDVT, WORK, LWORK, RESULT( 10 ) )
               END IF
               RESULT( 11 ) = ZERO
               DO 90 I = 1, MNMIN - 1
                  IF( SSAV( I ).LT.SSAV( I+1 ) ) RESULT( 11 ) = ULPINV                   IF( SSAV( I ).LT.ZERO ) RESULT( 11 ) = ULPINV
   90          CONTINUE
               IF( MNMIN.GE.1 ) THEN
                  IF( SSAV( MNMIN ).LT.ZERO ) RESULT( 11 ) = ULPINV
               END IF
*
*              Do partial SVDs, comparing to SSAV, USAV, and VTSAV
*
               RESULT( 12 ) = ZERO
               RESULT( 13 ) = ZERO
               RESULT( 14 ) = ZERO
               DO 110 IJQ = 0, 2
                  JOBQ = CJOB( IJQ+1 )
                  CALL DLACPY( 'F', M, N, ASAV, LDA, A, LDA )
                  SRNAMT = 'DGESDD'
                  CALL DGESDD( JOBQ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LSWORK, IWORK, IINFO )
*
*                 Compare U
*
                  DIF = ZERO
                  IF( M.GT.0 .AND. N.GT.0 ) THEN
                     IF( IJQ.EQ.1 ) THEN
                        IF( M.GE.N ) THEN
                           CALL DORT03( 'C', M, MNMIN, M, MNMIN, USAV, LDU, A, LDA, WORK, LWORK, DIF, INFO )
                        ELSE
                           CALL DORT03( 'C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, DIF, INFO )
                        END IF
                     ELSE IF( IJQ.EQ.2 ) THEN
                        CALL DORT03( 'C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, DIF, INFO )
                     END IF
                  END IF
                  RESULT( 12 ) = MAX( RESULT( 12 ), DIF )
*
*                 Compare VT
*
                  DIF = ZERO
                  IF( M.GT.0 .AND. N.GT.0 ) THEN
                     IF( IJQ.EQ.1 ) THEN
                        IF( M.GE.N ) THEN
                           CALL DORT03( 'R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, DIF, INFO )
                        ELSE
                           CALL DORT03( 'R', N, MNMIN, N, MNMIN, VTSAV, LDVT, A, LDA, WORK, LWORK, DIF, INFO )
                        END IF
                     ELSE IF( IJQ.EQ.2 ) THEN
                        CALL DORT03( 'R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, DIF, INFO )
                     END IF
                  END IF
                  RESULT( 13 ) = MAX( RESULT( 13 ), DIF )
*
*                 Compare S
*
                  DIF = ZERO
                  DIV = MAX( MNMIN*ULP*S( 1 ), UNFL )
                  DO 100 I = 1, MNMIN - 1
                     IF( SSAV( I ).LT.SSAV( I+1 ) ) DIF = ULPINV                      IF( SSAV( I ).LT.ZERO ) DIF = ULPINV
                     DIF = MAX( DIF, ABS( SSAV( I )-S( I ) ) / DIV )
  100             CONTINUE
                  RESULT( 14 ) = MAX( RESULT( 14 ), DIF )
  110          CONTINUE
*
*              Test DGESVDQ
*              Note: DGESVDQ only works for M >= N
*
               RESULT( 36 ) = ZERO
               RESULT( 37 ) = ZERO
               RESULT( 38 ) = ZERO
               RESULT( 39 ) = ZERO
*
               IF( M.GE.N ) THEN
                  IWTMP = 5*MNMIN*MNMIN + 9*MNMIN + MAX( M, N )
                  LSWORK = IWTMP + ( IWS-1 )*( LWORK-IWTMP ) / 3
                  LSWORK = MIN( LSWORK, LWORK )
                  LSWORK = MAX( LSWORK, 1 )
                  IF( IWS.EQ.4 ) LSWORK = LWORK
*
                  CALL DLACPY( 'F', M, N, ASAV, LDA, A, LDA )
                  SRNAMT = 'DGESVDQ'
*
                  LRWORK = 2
                  LIWORK = MAX( N, 1 )
                  CALL DGESVDQ( 'H', 'N', 'N', 'A', 'A',  M, N, A, LDA, SSAV, USAV, LDU, VTSAV, LDVT, NUMRANK, IWORK, LIWORK, WORK, LWORK, RWORK, LRWORK, IINFO )
*
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9995 )'DGESVDQ', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                     INFO = ABS( IINFO )
                     RETURN
                  END IF
*
*                 Do tests 36--39
*
                  CALL DBDT01( M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RESULT( 36 ) )
                  IF( M.NE.0 .AND. N.NE.0 ) THEN
                     CALL DORT01( 'Columns', M, M, USAV, LDU, WORK, LWORK, RESULT( 37 ) )                      CALL DORT01( 'Rows', N, N, VTSAV, LDVT, WORK, LWORK, RESULT( 38 ) )
                  END IF
                  RESULT( 39 ) = ZERO
                  DO 199 I = 1, MNMIN - 1
                     IF( SSAV( I ).LT.SSAV( I+1 ) ) RESULT( 39 ) = ULPINV                      IF( SSAV( I ).LT.ZERO ) RESULT( 39 ) = ULPINV
  199             CONTINUE
                  IF( MNMIN.GE.1 ) THEN
                     IF( SSAV( MNMIN ).LT.ZERO ) RESULT( 39 ) = ULPINV
                  END IF
               END IF
*
*              Test DGESVJ
*              Note: DGESVJ only works for M >= N
*
               RESULT( 15 ) = ZERO
               RESULT( 16 ) = ZERO
               RESULT( 17 ) = ZERO
               RESULT( 18 ) = ZERO
*
               IF( M.GE.N ) THEN
                  IWTMP = 5*MNMIN*MNMIN + 9*MNMIN + MAX( M, N )
                  LSWORK = IWTMP + ( IWS-1 )*( LWORK-IWTMP ) / 3
                  LSWORK = MIN( LSWORK, LWORK )
                  LSWORK = MAX( LSWORK, 1 )
                  IF( IWS.EQ.4 ) LSWORK = LWORK
*
                  CALL DLACPY( 'F', M, N, ASAV, LDA, USAV, LDA )
                  SRNAMT = 'DGESVJ'
                  CALL DGESVJ( 'G', 'U', 'V', M, N, USAV, LDA, SSAV, 0, A, LDVT, WORK, LWORK, INFO )
*
*                 DGESVJ returns V not VT
*
                  DO J=1,N
                     DO I=1,N
                        VTSAV(J,I) = A(I,J)
                     END DO
                  END DO
*
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9995 )'GESVJ', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                     INFO = ABS( IINFO )
                     RETURN
                  END IF
*
*                 Do tests 15--18
*
                  CALL DBDT01( M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RESULT( 15 ) )
                  IF( M.NE.0 .AND. N.NE.0 ) THEN
                     CALL DORT01( 'Columns', M, M, USAV, LDU, WORK, LWORK, RESULT( 16 ) )                      CALL DORT01( 'Rows', N, N, VTSAV, LDVT, WORK, LWORK, RESULT( 17 ) )
                  END IF
                  RESULT( 18 ) = ZERO
                  DO 120 I = 1, MNMIN - 1
                     IF( SSAV( I ).LT.SSAV( I+1 ) ) RESULT( 18 ) = ULPINV                      IF( SSAV( I ).LT.ZERO ) RESULT( 18 ) = ULPINV
  120             CONTINUE
                  IF( MNMIN.GE.1 ) THEN
                     IF( SSAV( MNMIN ).LT.ZERO ) RESULT( 18 ) = ULPINV
                  END IF
               END IF
*
*              Test DGEJSV
*              Note: DGEJSV only works for M >= N
*
               RESULT( 19 ) = ZERO
               RESULT( 20 ) = ZERO
               RESULT( 21 ) = ZERO
               RESULT( 22 ) = ZERO
               IF( M.GE.N ) THEN
                  IWTMP = 5*MNMIN*MNMIN + 9*MNMIN + MAX( M, N )
                  LSWORK = IWTMP + ( IWS-1 )*( LWORK-IWTMP ) / 3
                  LSWORK = MIN( LSWORK, LWORK )
                  LSWORK = MAX( LSWORK, 1 )
                  IF( IWS.EQ.4 ) LSWORK = LWORK
*
                  CALL DLACPY( 'F', M, N, ASAV, LDA, VTSAV, LDA )
                  SRNAMT = 'DGEJSV'
                  CALL DGEJSV( 'G', 'U', 'V', 'R', 'N', 'N', M, N, VTSAV, LDA, SSAV, USAV, LDU, A, LDVT, WORK, LWORK, IWORK, INFO )
*
*                 DGEJSV returns V not VT
*
                  DO 140 J=1,N
                     DO 130 I=1,N
                        VTSAV(J,I) = A(I,J)
  130                END DO
  140             END DO
*
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9995 )'GEJSV', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                     INFO = ABS( IINFO )
                     RETURN
                  END IF
*
*                 Do tests 19--22
*
                  CALL DBDT01( M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RESULT( 19 ) )
                  IF( M.NE.0 .AND. N.NE.0 ) THEN
                     CALL DORT01( 'Columns', M, M, USAV, LDU, WORK, LWORK, RESULT( 20 ) )                      CALL DORT01( 'Rows', N, N, VTSAV, LDVT, WORK, LWORK, RESULT( 21 ) )
                  END IF
                  RESULT( 22 ) = ZERO
                  DO 150 I = 1, MNMIN - 1
                     IF( SSAV( I ).LT.SSAV( I+1 ) ) RESULT( 22 ) = ULPINV                      IF( SSAV( I ).LT.ZERO ) RESULT( 22 ) = ULPINV
  150             CONTINUE
                  IF( MNMIN.GE.1 ) THEN
                     IF( SSAV( MNMIN ).LT.ZERO ) RESULT( 22 ) = ULPINV
                  END IF
               END IF
*
*              Test DGESVDX
*
               CALL DLACPY( 'F', M, N, ASAV, LDA, A, LDA )
               CALL DGESVDX( 'V', 'V', 'A', M, N, A, LDA, VL, VU, IL, IU, NS, SSAV, USAV, LDU, VTSAV, LDVT, WORK, LWORK, IWORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUT, FMT = 9995 )'GESVDX', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               END IF
*
*              Do tests 23--29
*
               RESULT( 23 ) = ZERO
               RESULT( 24 ) = ZERO
               RESULT( 25 ) = ZERO
               CALL DBDT01( M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RESULT( 23 ) )
               IF( M.NE.0 .AND. N.NE.0 ) THEN
                  CALL DORT01( 'Columns', M, M, USAV, LDU, WORK, LWORK, RESULT( 24 ) )                   CALL DORT01( 'Rows', N, N, VTSAV, LDVT, WORK, LWORK, RESULT( 25 ) )
               END IF
               RESULT( 26 ) = ZERO
               DO 160 I = 1, MNMIN - 1
                  IF( SSAV( I ).LT.SSAV( I+1 ) ) RESULT( 26 ) = ULPINV                   IF( SSAV( I ).LT.ZERO ) RESULT( 26 ) = ULPINV
  160          CONTINUE
               IF( MNMIN.GE.1 ) THEN
                  IF( SSAV( MNMIN ).LT.ZERO ) RESULT( 26 ) = ULPINV
               END IF
*
*              Do partial SVDs, comparing to SSAV, USAV, and VTSAV
*
               RESULT( 27 ) = ZERO
               RESULT( 28 ) = ZERO
               RESULT( 29 ) = ZERO
               DO 180 IJU = 0, 1
                  DO 170 IJVT = 0, 1
                     IF( ( IJU.EQ.0 .AND. IJVT.EQ.0 ) .OR. ( IJU.EQ.1 .AND. IJVT.EQ.1 ) )GO TO 170
                     JOBU = CJOBV( IJU+1 )
                     JOBVT = CJOBV( IJVT+1 )
                     RANGE = CJOBR( 1 )
                     CALL DLACPY( 'F', M, N, ASAV, LDA, A, LDA )
                     CALL DGESVDX( JOBU, JOBVT, RANGE, M, N, A, LDA, VL, VU, IL, IU, NS, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, IINFO )
*
*                    Compare U
*
                     DIF = ZERO
                     IF( M.GT.0 .AND. N.GT.0 ) THEN
                        IF( IJU.EQ.1 ) THEN
                           CALL DORT03( 'C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, DIF, IINFO )
                        END IF
                     END IF
                     RESULT( 27 ) = MAX( RESULT( 27 ), DIF )
*
*                    Compare VT
*
                     DIF = ZERO
                     IF( M.GT.0 .AND. N.GT.0 ) THEN
                        IF( IJVT.EQ.1 ) THEN
                           CALL DORT03( 'R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, DIF, IINFO )
                        END IF
                     END IF
                     RESULT( 28 ) = MAX( RESULT( 28 ), DIF )
*
*                    Compare S
*
                     DIF = ZERO
                     DIV = MAX( MNMIN*ULP*S( 1 ), UNFL )
                     DO 190 I = 1, MNMIN - 1
                        IF( SSAV( I ).LT.SSAV( I+1 ) ) DIF = ULPINV                         IF( SSAV( I ).LT.ZERO ) DIF = ULPINV
                        DIF = MAX( DIF, ABS( SSAV( I )-S( I ) ) / DIV )
  190                CONTINUE
                     RESULT( 29 ) = MAX( RESULT( 29 ), DIF )
  170             CONTINUE
  180          CONTINUE
*
*              Do tests 30--32: DGESVDX( 'V', 'V', 'I' )
*
               DO 200 I = 1, 4
                  ISEED2( I ) = ISEED( I )
  200          CONTINUE
               IF( MNMIN.LE.1 ) THEN
                  IL = 1
                  IU = MAX( 1, MNMIN )
               ELSE
                  IL = 1 + INT( ( MNMIN-1 )*DLARND( 1, ISEED2 ) )
                  IU = 1 + INT( ( MNMIN-1 )*DLARND( 1, ISEED2 ) )
                  IF( IU.LT.IL ) THEN
                     ITEMP = IU
                     IU = IL
                     IL = ITEMP
                  END IF
               END IF
               CALL DLACPY( 'F', M, N, ASAV, LDA, A, LDA )
               CALL DGESVDX( 'V', 'V', 'I', M, N, A, LDA, VL, VU, IL, IU, NSI, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUT, FMT = 9995 )'GESVDX', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               END IF
*
               RESULT( 30 ) = ZERO
               RESULT( 31 ) = ZERO
               RESULT( 32 ) = ZERO
               CALL DBDT05( M, N, ASAV, LDA, S, NSI, U, LDU, VT, LDVT, WORK, RESULT( 30 ) )                CALL DORT01( 'Columns', M, NSI, U, LDU, WORK, LWORK, RESULT( 31 ) )                CALL DORT01( 'Rows', NSI, N, VT, LDVT, WORK, LWORK, RESULT( 32 ) )
*
*              Do tests 33--35: DGESVDX( 'V', 'V', 'V' )
*
               IF( MNMIN.GT.0 .AND. NSI.GT.1 ) THEN
                  IF( IL.NE.1 ) THEN
                     VU = SSAV( IL ) + MAX( HALF*ABS( SSAV( IL )-SSAV( IL-1 ) ), ULP*ANORM, TWO*RTUNFL )
                  ELSE
                     VU = SSAV( 1 ) + MAX( HALF*ABS( SSAV( NS )-SSAV( 1 ) ), ULP*ANORM, TWO*RTUNFL )
                  END IF
                  IF( IU.NE.NS ) THEN
                     VL = SSAV( IU ) - MAX( ULP*ANORM, TWO*RTUNFL, HALF*ABS( SSAV( IU+1 )-SSAV( IU ) ) )
                  ELSE
                     VL = SSAV( NS ) - MAX( ULP*ANORM, TWO*RTUNFL, HALF*ABS( SSAV( NS )-SSAV( 1 ) ) )
                  END IF
                  VL = MAX( VL,ZERO )
                  VU = MAX( VU,ZERO )
                  IF( VL.GE.VU ) VU = MAX( VU*2, VU+VL+HALF )
               ELSE
                  VL = ZERO
                  VU = ONE
               END IF
               CALL DLACPY( 'F', M, N, ASAV, LDA, A, LDA )
               CALL DGESVDX( 'V', 'V', 'V', M, N, A, LDA, VL, VU, IL, IU, NSV, S, U, LDU, VT, LDVT, WORK, LWORK, IWORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUT, FMT = 9995 )'GESVDX', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               END IF
*
               RESULT( 33 ) = ZERO
               RESULT( 34 ) = ZERO
               RESULT( 35 ) = ZERO
               CALL DBDT05( M, N, ASAV, LDA, S, NSV, U, LDU, VT, LDVT, WORK, RESULT( 33 ) )                CALL DORT01( 'Columns', M, NSV, U, LDU, WORK, LWORK, RESULT( 34 ) )                CALL DORT01( 'Rows', NSV, N, VT, LDVT, WORK, LWORK, RESULT( 35 ) )
*
*              End of Loop -- Check for RESULT(j) > THRESH
*
               DO 210 J = 1, 39
                  IF( RESULT( J ).GE.THRESH ) THEN
                     IF( NFAIL.EQ.0 ) THEN
                        WRITE( NOUT, FMT = 9999 )
                        WRITE( NOUT, FMT = 9998 )
                     END IF
                     WRITE( NOUT, FMT = 9997 )M, N, JTYPE, IWS, IOLDSD, J, RESULT( J )
                     NFAIL = NFAIL + 1
                  END IF
  210          CONTINUE
               NTEST = NTEST + 39
  220       CONTINUE
  230    CONTINUE
  240 CONTINUE
*
*     Summary
*
      CALL ALASVM( PATH, NOUT, NFAIL, NTEST, 0 )
*
 9999 FORMAT( ' SVD -- Real Singular Value Decomposition Driver ',
     $      / ' Matrix types (see DDRVBD for details):',
     $      / / ' 1 = Zero matrix', / ' 2 = Identity matrix',
     $      / ' 3 = Evenly spaced singular values near 1',
     $      / ' 4 = Evenly spaced singular values near underflow',
     $      / ' 5 = Evenly spaced singular values near overflow', / /
     $      ' Tests performed: ( A is dense, U and V are orthogonal,',
     $      / 19X, ' S is an array, and Upartial, VTpartial, and',
     $      / 19X, ' Spartial are partially computed U, VT and S),', / )
 9998 FORMAT( ' 1 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ',
     $      / ' 2 = | I - U**T U | / ( M ulp ) ',
     $      / ' 3 = | I - VT VT**T | / ( N ulp ) ',
     $      / ' 4 = 0 if S contains min(M,N) nonnegative values in',
     $      ' decreasing order, else 1/ulp',
     $      / ' 5 = | U - Upartial | / ( M ulp )',
     $      / ' 6 = | VT - VTpartial | / ( N ulp )',
     $      / ' 7 = | S - Spartial | / ( min(M,N) ulp |S| )',
     $      / ' 8 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ',
     $      / ' 9 = | I - U**T U | / ( M ulp ) ',
     $      / '10 = | I - VT VT**T | / ( N ulp ) ',
     $      / '11 = 0 if S contains min(M,N) nonnegative values in',
     $      ' decreasing order, else 1/ulp',
     $      / '12 = | U - Upartial | / ( M ulp )',
     $      / '13 = | VT - VTpartial | / ( N ulp )',
     $      / '14 = | S - Spartial | / ( min(M,N) ulp |S| )',
     $      / '15 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ',
     $      / '16 = | I - U**T U | / ( M ulp ) ',
     $      / '17 = | I - VT VT**T | / ( N ulp ) ',
     $      / '18 = 0 if S contains min(M,N) nonnegative values in',
     $      ' decreasing order, else 1/ulp',
     $      / '19 = | U - Upartial | / ( M ulp )',
     $      / '20 = | VT - VTpartial | / ( N ulp )',
     $      / '21 = | S - Spartial | / ( min(M,N) ulp |S| )',
     $      / '22 = 0 if S contains min(M,N) nonnegative values in',
     $      ' decreasing order, else 1/ulp',
     $      / '23 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ),',
     $      ' DGESVDX(V,V,A) ',
     $      / '24 = | I - U**T U | / ( M ulp ) ',
     $      / '25 = | I - VT VT**T | / ( N ulp ) ',
     $      / '26 = 0 if S contains min(M,N) nonnegative values in',
     $      ' decreasing order, else 1/ulp',
     $      / '27 = | U - Upartial | / ( M ulp )',
     $      / '28 = | VT - VTpartial | / ( N ulp )',
     $      / '29 = | S - Spartial | / ( min(M,N) ulp |S| )',
     $      / '30 = | U**T A VT**T - diag(S) | / ( |A| max(M,N) ulp ),',
     $      ' DGESVDX(V,V,I) ',
     $      / '31 = | I - U**T U | / ( M ulp ) ',
     $      / '32 = | I - VT VT**T | / ( N ulp ) ',
     $      / '33 = | U**T A VT**T - diag(S) | / ( |A| max(M,N) ulp ),',
     $      ' DGESVDX(V,V,V) ',
     $      / '34 = | I - U**T U | / ( M ulp ) ',
     $      / '35 = | I - VT VT**T | / ( N ulp ) ',
     $      ' DGESVDQ(H,N,N,A,A',
     $      / '36 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ',
     $      / '37 = | I - U**T U | / ( M ulp ) ',
     $      / '38 = | I - VT VT**T | / ( N ulp ) ',
     $      / '39 = 0 if S contains min(M,N) nonnegative values in',
     $      ' decreasing order, else 1/ulp',
     $      / / )
 9997 FORMAT( ' M=', I5, ', N=', I5, ', type ', I1, ', IWS=', I1,
     $      ', seed=', 4( I4, ',' ), ' test(', I2, ')=', G11.4 )
 9996 FORMAT( ' DDRVBD: ', A, ' returned INFO=', I6, '.', / 9X, 'M=',
     $      I6, ', N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ),
     $      I5, ')' )
 9995 FORMAT( ' DDRVBD: ', A, ' returned INFO=', I6, '.', / 9X, 'M=',
     $      I6, ', N=', I6, ', JTYPE=', I6, ', LSWORK=', I6, / 9X,
     $      'ISEED=(', 3( I5, ',' ), I5, ')' )
*
      RETURN
*
*     End of DDRVBD
*
      END
