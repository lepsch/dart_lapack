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

      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZDRVBD', -INFO )
         RETURN
      END IF

      // Quick return if nothing to do

      IF( NSIZES.EQ.0 .OR. NTYPES.EQ.0 ) RETURN

      // More Important constants

      UNFL = DLAMCH( 'S' )
      OVFL = ONE / UNFL
      ULP = DLAMCH( 'E' )
      ULPINV = ONE / ULP
      RTUNFL = SQRT( UNFL )

      // Loop over sizes, types

      NERRS = 0

      DO 230 JSIZE = 1, NSIZES
         M = MM( JSIZE )
         N = NN( JSIZE )
         MNMIN = MIN( M, N )

         IF( NSIZES.NE.1 ) THEN
            MTYPES = MIN( MAXTYP, NTYPES )
         ELSE
            MTYPES = MIN( MAXTYP+1, NTYPES )
         END IF

         DO 220 JTYPE = 1, MTYPES
            IF( .NOT.DOTYPE( JTYPE ) ) GO TO 220
            NTEST = 0

            DO 20 J = 1, 4
               IOLDSD( J ) = ISEED( J )
   20       CONTINUE

            // Compute "A"

            IF( MTYPES.GT.MAXTYP ) GO TO 50

            IF( JTYPE.EQ.1 ) THEN

               // Zero matrix

               CALL ZLASET( 'Full', M, N, CZERO, CZERO, A, LDA )
               DO 30 I = 1, MIN( M, N )
                  S( I ) = ZERO
   30          CONTINUE

            ELSE IF( JTYPE.EQ.2 ) THEN

               // Identity matrix

               CALL ZLASET( 'Full', M, N, CZERO, CONE, A, LDA )
               DO 40 I = 1, MIN( M, N )
                  S( I ) = ONE
   40          CONTINUE

            ELSE

               // (Scaled) random matrix

               IF( JTYPE.EQ.3 ) ANORM = ONE                IF( JTYPE.EQ.4 ) ANORM = UNFL / ULP                IF( JTYPE.EQ.5 ) ANORM = OVFL*ULP                CALL ZLATMS( M, N, 'U', ISEED, 'N', S, 4, DBLE( MNMIN ), ANORM, M-1, N-1, 'N', A, LDA, WORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9996 )'Generator', IINFO, M, N, JTYPE, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               END IF
            END IF

   50       CONTINUE
            CALL ZLACPY( 'F', M, N, A, LDA, ASAV, LDA )

            // Do for minimal and adequate (for blocking) workspace

            DO 210 IWSPC = 1, 4

               // Test for ZGESVD

               IWTMP = 2*MIN( M, N )+MAX( M, N )
               LSWORK = IWTMP + ( IWSPC-1 )*( LWORK-IWTMP ) / 3
               LSWORK = MIN( LSWORK, LWORK )
               LSWORK = MAX( LSWORK, 1 )
               IF( IWSPC.EQ.4 ) LSWORK = LWORK

               DO 60 J = 1, 35
                  RESULT( J ) = -ONE
   60          CONTINUE

               // Factorize A

               IF( IWSPC.GT.1 ) CALL ZLACPY( 'F', M, N, ASAV, LDA, A, LDA )
               SRNAMT = 'ZGESVD'
               CALL ZGESVD( 'A', 'A', M, N, A, LDA, SSAV, USAV, LDU, VTSAV, LDVT, WORK, LSWORK, RWORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9995 )'GESVD', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               END IF

               // Do tests 1--4

               CALL ZBDT01( M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RWORK, RESULT( 1 ) )
               IF( M.NE.0 .AND. N.NE.0 ) THEN
                  CALL ZUNT01( 'Columns', MNMIN, M, USAV, LDU, WORK, LWORK, RWORK, RESULT( 2 ) )                   CALL ZUNT01( 'Rows', MNMIN, N, VTSAV, LDVT, WORK, LWORK, RWORK, RESULT( 3 ) )
               END IF
               RESULT( 4 ) = 0
               DO 70 I = 1, MNMIN - 1
                  IF( SSAV( I ).LT.SSAV( I+1 ) ) RESULT( 4 ) = ULPINV                   IF( SSAV( I ).LT.ZERO ) RESULT( 4 ) = ULPINV
   70          CONTINUE
               IF( MNMIN.GE.1 ) THEN
                  IF( SSAV( MNMIN ).LT.ZERO ) RESULT( 4 ) = ULPINV
               END IF

               // Do partial SVDs, comparing to SSAV, USAV, and VTSAV

               RESULT( 5 ) = ZERO
               RESULT( 6 ) = ZERO
               RESULT( 7 ) = ZERO
               DO 100 IJU = 0, 3
                  DO 90 IJVT = 0, 3
                     IF( ( IJU.EQ.3 .AND. IJVT.EQ.3 ) .OR. ( IJU.EQ.1 .AND. IJVT.EQ.1 ) )GO TO 90
                     JOBU = CJOB( IJU+1 )
                     JOBVT = CJOB( IJVT+1 )
                     CALL ZLACPY( 'F', M, N, ASAV, LDA, A, LDA )
                     SRNAMT = 'ZGESVD'
                     CALL ZGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LSWORK, RWORK, IINFO )

                     // Compare U

                     DIF = ZERO
                     IF( M.GT.0 .AND. N.GT.0 ) THEN
                        IF( IJU.EQ.1 ) THEN
                           CALL ZUNT03( 'C', M, MNMIN, M, MNMIN, USAV, LDU, A, LDA, WORK, LWORK, RWORK, DIF, IINFO )
                        ELSE IF( IJU.EQ.2 ) THEN
                           CALL ZUNT03( 'C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, RWORK, DIF, IINFO )
                        ELSE IF( IJU.EQ.3 ) THEN
                           CALL ZUNT03( 'C', M, M, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, RWORK, DIF, IINFO )
                        END IF
                     END IF
                     RESULT( 5 ) = MAX( RESULT( 5 ), DIF )

                     // Compare VT

                     DIF = ZERO
                     IF( M.GT.0 .AND. N.GT.0 ) THEN
                        IF( IJVT.EQ.1 ) THEN
                           CALL ZUNT03( 'R', N, MNMIN, N, MNMIN, VTSAV, LDVT, A, LDA, WORK, LWORK, RWORK, DIF, IINFO )
                        ELSE IF( IJVT.EQ.2 ) THEN
                           CALL ZUNT03( 'R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, RWORK, DIF, IINFO )
                        ELSE IF( IJVT.EQ.3 ) THEN
                           CALL ZUNT03( 'R', N, N, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, RWORK, DIF, IINFO )
                        END IF
                     END IF
                     RESULT( 6 ) = MAX( RESULT( 6 ), DIF )

                     // Compare S

                     DIF = ZERO
                     DIV = MAX( DBLE( MNMIN )*ULP*S( 1 ), DLAMCH( 'Safe minimum' ) )
                     DO 80 I = 1, MNMIN - 1
                        IF( SSAV( I ).LT.SSAV( I+1 ) ) DIF = ULPINV                         IF( SSAV( I ).LT.ZERO ) DIF = ULPINV
                        DIF = MAX( DIF, ABS( SSAV( I )-S( I ) ) / DIV )
   80                CONTINUE
                     RESULT( 7 ) = MAX( RESULT( 7 ), DIF )
   90             CONTINUE
  100          CONTINUE

               // Test for ZGESDD

               IWTMP = 2*MNMIN*MNMIN + 2*MNMIN + MAX( M, N )
               LSWORK = IWTMP + ( IWSPC-1 )*( LWORK-IWTMP ) / 3
               LSWORK = MIN( LSWORK, LWORK )
               LSWORK = MAX( LSWORK, 1 )
               IF( IWSPC.EQ.4 ) LSWORK = LWORK

               // Factorize A

               CALL ZLACPY( 'F', M, N, ASAV, LDA, A, LDA )
               SRNAMT = 'ZGESDD'
               CALL ZGESDD( 'A', M, N, A, LDA, SSAV, USAV, LDU, VTSAV, LDVT, WORK, LSWORK, RWORK, IWORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9995 )'GESDD', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               END IF

               // Do tests 1--4

               CALL ZBDT01( M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RWORK, RESULT( 8 ) )
               IF( M.NE.0 .AND. N.NE.0 ) THEN
                  CALL ZUNT01( 'Columns', MNMIN, M, USAV, LDU, WORK, LWORK, RWORK, RESULT( 9 ) )                   CALL ZUNT01( 'Rows', MNMIN, N, VTSAV, LDVT, WORK, LWORK, RWORK, RESULT( 10 ) )
               END IF
               RESULT( 11 ) = 0
               DO 110 I = 1, MNMIN - 1
                  IF( SSAV( I ).LT.SSAV( I+1 ) ) RESULT( 11 ) = ULPINV                   IF( SSAV( I ).LT.ZERO ) RESULT( 11 ) = ULPINV
  110          CONTINUE
               IF( MNMIN.GE.1 ) THEN
                  IF( SSAV( MNMIN ).LT.ZERO ) RESULT( 11 ) = ULPINV
               END IF

               // Do partial SVDs, comparing to SSAV, USAV, and VTSAV

               RESULT( 12 ) = ZERO
               RESULT( 13 ) = ZERO
               RESULT( 14 ) = ZERO
               DO 130 IJQ = 0, 2
                  JOBQ = CJOB( IJQ+1 )
                  CALL ZLACPY( 'F', M, N, ASAV, LDA, A, LDA )
                  SRNAMT = 'ZGESDD'
                  CALL ZGESDD( JOBQ, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LSWORK, RWORK, IWORK, IINFO )

                  // Compare U

                  DIF = ZERO
                  IF( M.GT.0 .AND. N.GT.0 ) THEN
                     IF( IJQ.EQ.1 ) THEN
                        IF( M.GE.N ) THEN
                           CALL ZUNT03( 'C', M, MNMIN, M, MNMIN, USAV, LDU, A, LDA, WORK, LWORK, RWORK, DIF, IINFO )
                        ELSE
                           CALL ZUNT03( 'C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, RWORK, DIF, IINFO )
                        END IF
                     ELSE IF( IJQ.EQ.2 ) THEN
                        CALL ZUNT03( 'C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, RWORK, DIF, IINFO )
                     END IF
                  END IF
                  RESULT( 12 ) = MAX( RESULT( 12 ), DIF )

                  // Compare VT

                  DIF = ZERO
                  IF( M.GT.0 .AND. N.GT.0 ) THEN
                     IF( IJQ.EQ.1 ) THEN
                        IF( M.GE.N ) THEN
                           CALL ZUNT03( 'R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, RWORK, DIF, IINFO )
                        ELSE
                           CALL ZUNT03( 'R', N, MNMIN, N, MNMIN, VTSAV, LDVT, A, LDA, WORK, LWORK, RWORK, DIF, IINFO )
                        END IF
                     ELSE IF( IJQ.EQ.2 ) THEN
                        CALL ZUNT03( 'R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, RWORK, DIF, IINFO )
                     END IF
                  END IF
                  RESULT( 13 ) = MAX( RESULT( 13 ), DIF )

                  // Compare S

                  DIF = ZERO
                  DIV = MAX( DBLE( MNMIN )*ULP*S( 1 ), DLAMCH( 'Safe minimum' ) )
                  DO 120 I = 1, MNMIN - 1
                     IF( SSAV( I ).LT.SSAV( I+1 ) ) DIF = ULPINV                      IF( SSAV( I ).LT.ZERO ) DIF = ULPINV
                     DIF = MAX( DIF, ABS( SSAV( I )-S( I ) ) / DIV )
  120             CONTINUE
                  RESULT( 14 ) = MAX( RESULT( 14 ), DIF )
  130          CONTINUE

               // Test ZGESVDQ
               // Note: ZGESVDQ only works for M >= N

               RESULT( 36 ) = ZERO
               RESULT( 37 ) = ZERO
               RESULT( 38 ) = ZERO
               RESULT( 39 ) = ZERO

               IF( M.GE.N ) THEN
                  IWTMP = 2*MNMIN*MNMIN + 2*MNMIN + MAX( M, N )
                  LSWORK = IWTMP + ( IWSPC-1 )*( LWORK-IWTMP ) / 3
                  LSWORK = MIN( LSWORK, LWORK )
                  LSWORK = MAX( LSWORK, 1 )
                  IF( IWSPC.EQ.4 ) LSWORK = LWORK

                  CALL ZLACPY( 'F', M, N, ASAV, LDA, A, LDA )
                  SRNAMT = 'ZGESVDQ'

                  LRWORK = MAX(2, M, 5*N)
                  LIWORK = MAX( N, 1 )
                  CALL ZGESVDQ( 'H', 'N', 'N', 'A', 'A',  M, N, A, LDA, SSAV, USAV, LDU, VTSAV, LDVT, NUMRANK, IWORK, LIWORK, WORK, LWORK, RWORK, LRWORK, IINFO )

                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUNIT, FMT = 9995 )'ZGESVDQ', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                     INFO = ABS( IINFO )
                     RETURN
                  END IF

                  // Do tests 36--39

                  CALL ZBDT01( M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RWORK, RESULT( 36 ) )
                  IF( M.NE.0 .AND. N.NE.0 ) THEN
                     CALL ZUNT01( 'Columns', M, M, USAV, LDU, WORK, LWORK, RWORK, RESULT( 37 ) )                      CALL ZUNT01( 'Rows', N, N, VTSAV, LDVT, WORK, LWORK, RWORK, RESULT( 38 ) )
                  END IF
                  RESULT( 39 ) = ZERO
                  DO 199 I = 1, MNMIN - 1
                     IF( SSAV( I ).LT.SSAV( I+1 ) ) RESULT( 39 ) = ULPINV                      IF( SSAV( I ).LT.ZERO ) RESULT( 39 ) = ULPINV
  199             CONTINUE
                  IF( MNMIN.GE.1 ) THEN
                     IF( SSAV( MNMIN ).LT.ZERO ) RESULT( 39 ) = ULPINV
                  END IF
               END IF

               // Test ZGESVJ
               // Note: ZGESVJ only works for M >= N

               RESULT( 15 ) = ZERO
               RESULT( 16 ) = ZERO
               RESULT( 17 ) = ZERO
               RESULT( 18 ) = ZERO

               IF( M.GE.N ) THEN
                  IWTMP = 2*MNMIN*MNMIN + 2*MNMIN + MAX( M, N )
                  LSWORK = IWTMP + ( IWSPC-1 )*( LWORK-IWTMP ) / 3
                  LSWORK = MIN( LSWORK, LWORK )
                  LSWORK = MAX( LSWORK, 1 )
                  LRWORK = MAX(6,N)
                  IF( IWSPC.EQ.4 ) LSWORK = LWORK

                  CALL ZLACPY( 'F', M, N, ASAV, LDA, USAV, LDA )
                  SRNAMT = 'ZGESVJ'
                  CALL ZGESVJ( 'G', 'U', 'V', M, N, USAV, LDA, SSAV, 0, A, LDVT, WORK, LWORK, RWORK, LRWORK, IINFO )

                  // ZGESVJ returns V not VH

                  DO J=1,N
                     DO I=1,N
                        VTSAV(J,I) = CONJG (A(I,J))
                     END DO
                  END DO

                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUNIT, FMT = 9995 )'GESVJ', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                     INFO = ABS( IINFO )
                     RETURN
                  END IF

                  // Do tests 15--18

                  CALL ZBDT01( M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RWORK, RESULT( 15 ) )
                  IF( M.NE.0 .AND. N.NE.0 ) THEN
                     CALL ZUNT01( 'Columns', M, M, USAV, LDU, WORK, LWORK, RWORK, RESULT( 16 ) )                      CALL ZUNT01( 'Rows', N, N, VTSAV, LDVT, WORK, LWORK, RWORK, RESULT( 17 ) )
                  END IF
                  RESULT( 18 ) = ZERO
                  DO 131 I = 1, MNMIN - 1
                     IF( SSAV( I ).LT.SSAV( I+1 ) ) RESULT( 18 ) = ULPINV                      IF( SSAV( I ).LT.ZERO ) RESULT( 18 ) = ULPINV
  131             CONTINUE
                  IF( MNMIN.GE.1 ) THEN
                     IF( SSAV( MNMIN ).LT.ZERO ) RESULT( 18 ) = ULPINV
                  END IF
               END IF

               // Test ZGEJSV
               // Note: ZGEJSV only works for M >= N

               RESULT( 19 ) = ZERO
               RESULT( 20 ) = ZERO
               RESULT( 21 ) = ZERO
               RESULT( 22 ) = ZERO
               IF( M.GE.N ) THEN
                  IWTMP = 2*MNMIN*MNMIN + 2*MNMIN + MAX( M, N )
                  LSWORK = IWTMP + ( IWSPC-1 )*( LWORK-IWTMP ) / 3
                  LSWORK = MIN( LSWORK, LWORK )
                  LSWORK = MAX( LSWORK, 1 )
                  IF( IWSPC.EQ.4 ) LSWORK = LWORK
                  LRWORK = MAX( 7, N + 2*M)

                  CALL ZLACPY( 'F', M, N, ASAV, LDA, VTSAV, LDA )
                  SRNAMT = 'ZGEJSV'
                  CALL ZGEJSV( 'G', 'U', 'V', 'R', 'N', 'N', M, N, VTSAV, LDA, SSAV, USAV, LDU, A, LDVT, WORK, LWORK, RWORK, LRWORK, IWORK, IINFO )

                  // ZGEJSV returns V not VH

                  DO 133 J=1,N
                     DO 132 I=1,N
                        VTSAV(J,I) = CONJG (A(I,J))
  132                END DO
  133             END DO

                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUNIT, FMT = 9995 )'GEJSV', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                     INFO = ABS( IINFO )
                     RETURN
                  END IF

                  // Do tests 19--22

                  CALL ZBDT01( M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RWORK, RESULT( 19 ) )
                  IF( M.NE.0 .AND. N.NE.0 ) THEN
                     CALL ZUNT01( 'Columns', M, M, USAV, LDU, WORK, LWORK, RWORK, RESULT( 20 ) )                      CALL ZUNT01( 'Rows', N, N, VTSAV, LDVT, WORK, LWORK, RWORK, RESULT( 21 ) )
                  END IF
                  RESULT( 22 ) = ZERO
                  DO 134 I = 1, MNMIN - 1
                     IF( SSAV( I ).LT.SSAV( I+1 ) ) RESULT( 22 ) = ULPINV                      IF( SSAV( I ).LT.ZERO ) RESULT( 22 ) = ULPINV
  134             CONTINUE
                  IF( MNMIN.GE.1 ) THEN
                     IF( SSAV( MNMIN ).LT.ZERO ) RESULT( 22 ) = ULPINV
                  END IF
               END IF

               // Test ZGESVDX

               // Factorize A

               CALL ZLACPY( 'F', M, N, ASAV, LDA, A, LDA )
               SRNAMT = 'ZGESVDX'
               CALL ZGESVDX( 'V', 'V', 'A', M, N, A, LDA, VL, VU, IL, IU, NS, SSAV, USAV, LDU, VTSAV, LDVT, WORK, LWORK, RWORK, IWORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9995 )'GESVDX', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               END IF

               // Do tests 1--4

               RESULT( 23 ) = ZERO
               RESULT( 24 ) = ZERO
               RESULT( 25 ) = ZERO
               CALL ZBDT01( M, N, 0, ASAV, LDA, USAV, LDU, SSAV, E, VTSAV, LDVT, WORK, RWORK, RESULT( 23 ) )
               IF( M.NE.0 .AND. N.NE.0 ) THEN
                  CALL ZUNT01( 'Columns', MNMIN, M, USAV, LDU, WORK, LWORK, RWORK, RESULT( 24 ) )                   CALL ZUNT01( 'Rows', MNMIN, N, VTSAV, LDVT, WORK, LWORK, RWORK, RESULT( 25 ) )
               END IF
               RESULT( 26 ) = ZERO
               DO 140 I = 1, MNMIN - 1
                  IF( SSAV( I ).LT.SSAV( I+1 ) ) RESULT( 26 ) = ULPINV                   IF( SSAV( I ).LT.ZERO ) RESULT( 26 ) = ULPINV
  140          CONTINUE
               IF( MNMIN.GE.1 ) THEN
                  IF( SSAV( MNMIN ).LT.ZERO ) RESULT( 26 ) = ULPINV
               END IF

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
                     CALL ZLACPY( 'F', M, N, ASAV, LDA, A, LDA )
                     SRNAMT = 'ZGESVDX'
                     CALL ZGESVDX( JOBU, JOBVT, 'A', M, N, A, LDA, VL, VU, IL, IU, NS, SSAV, U, LDU, VT, LDVT, WORK, LWORK, RWORK, IWORK, IINFO )

                     // Compare U

                     DIF = ZERO
                     IF( M.GT.0 .AND. N.GT.0 ) THEN
                        IF( IJU.EQ.1 ) THEN
                           CALL ZUNT03( 'C', M, MNMIN, M, MNMIN, USAV, LDU, U, LDU, WORK, LWORK, RWORK, DIF, IINFO )
                        END IF
                     END IF
                     RESULT( 27 ) = MAX( RESULT( 27 ), DIF )

                     // Compare VT

                     DIF = ZERO
                     IF( M.GT.0 .AND. N.GT.0 ) THEN
                        IF( IJVT.EQ.1 ) THEN
                           CALL ZUNT03( 'R', N, MNMIN, N, MNMIN, VTSAV, LDVT, VT, LDVT, WORK, LWORK, RWORK, DIF, IINFO )
                        END IF
                     END IF
                     RESULT( 28 ) = MAX( RESULT( 28 ), DIF )

                     // Compare S

                     DIF = ZERO
                     DIV = MAX( DBLE( MNMIN )*ULP*S( 1 ), DLAMCH( 'Safe minimum' ) )
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
               CALL ZLACPY( 'F', M, N, ASAV, LDA, A, LDA )
               SRNAMT = 'ZGESVDX'
               CALL ZGESVDX( 'V', 'V', 'I', M, N, A, LDA, VL, VU, IL, IU, NSI, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, IWORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9995 )'GESVDX', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               END IF

               RESULT( 30 ) = ZERO
               RESULT( 31 ) = ZERO
               RESULT( 32 ) = ZERO
               CALL ZBDT05( M, N, ASAV, LDA, S, NSI, U, LDU, VT, LDVT, WORK, RESULT( 30 ) )
               IF( M.NE.0 .AND. N.NE.0 ) THEN
                  CALL ZUNT01( 'Columns', M, NSI, U, LDU, WORK, LWORK, RWORK, RESULT( 31 ) )                   CALL ZUNT01( 'Rows', NSI, N, VT, LDVT, WORK, LWORK, RWORK, RESULT( 32 ) )
               END IF

               // Do tests 11--13

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
               CALL ZLACPY( 'F', M, N, ASAV, LDA, A, LDA )
               SRNAMT = 'ZGESVDX'
               CALL ZGESVDX( 'V', 'V', 'V', M, N, A, LDA, VL, VU, IL, IU, NSV, S, U, LDU, VT, LDVT, WORK, LWORK, RWORK, IWORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUNIT, FMT = 9995 )'GESVDX', IINFO, M, N, JTYPE, LSWORK, IOLDSD
                  INFO = ABS( IINFO )
                  RETURN
               END IF

               RESULT( 33 ) = ZERO
               RESULT( 34 ) = ZERO
               RESULT( 35 ) = ZERO
               CALL ZBDT05( M, N, ASAV, LDA, S, NSV, U, LDU, VT, LDVT, WORK, RESULT( 33 ) )
               IF( M.NE.0 .AND. N.NE.0 ) THEN
                  CALL ZUNT01( 'Columns', M, NSV, U, LDU, WORK, LWORK, RWORK, RESULT( 34 ) )                   CALL ZUNT01( 'Rows', NSV, N, VT, LDVT, WORK, LWORK, RWORK, RESULT( 35 ) )
               END IF

               // End of Loop -- Check for RESULT(j) > THRESH

               NTEST = 0
               NFAIL = 0
               DO 190 J = 1, 39
                  IF( RESULT( J ).GE.ZERO ) NTEST = NTEST + 1                   IF( RESULT( J ).GE.THRESH ) NFAIL = NFAIL + 1
  190          CONTINUE

               IF( NFAIL.GT.0 ) NTESTF = NTESTF + 1
               IF( NTESTF.EQ.1 ) THEN
                  WRITE( NOUNIT, FMT = 9999 )
                  WRITE( NOUNIT, FMT = 9998 )THRESH
                  NTESTF = 2
               END IF

               DO 200 J = 1, 39
                  IF( RESULT( J ).GE.THRESH ) THEN
                     WRITE( NOUNIT, FMT = 9997 )M, N, JTYPE, IWSPC, IOLDSD, J, RESULT( J )
                  END IF
  200          CONTINUE

               NERRS = NERRS + NFAIL
               NTESTT = NTESTT + NTEST

  210       CONTINUE

  220    CONTINUE
  230 CONTINUE

      // Summary

      CALL ALASVM( 'ZBD', NOUNIT, NERRS, NTESTT, 0 )

 9999 FORMAT( ' SVD -- Complex Singular Value Decomposition Driver ',
     $      / ' Matrix types (see ZDRVBD for details):',
     $      / / ' 1 = Zero matrix', / ' 2 = Identity matrix',
     $      / ' 3 = Evenly spaced singular values near 1',
     $      / ' 4 = Evenly spaced singular values near underflow',
     $      / ' 5 = Evenly spaced singular values near overflow',
     $      / / ' Tests performed: ( A is dense, U and V are unitary,',
     $      / 19X, ' S is an array, and Upartial, VTpartial, and',
     $      / 19X, ' Spartial are partially computed U, VT and S),', / )
 9998 FORMAT( ' Tests performed with Test Threshold = ', F8.2,
     $      / ' ZGESVD: ', /
     $      ' 1 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ',
     $      / ' 2 = | I - U**T U | / ( M ulp ) ',
     $      / ' 3 = | I - VT VT**T | / ( N ulp ) ',
     $      / ' 4 = 0 if S contains min(M,N) nonnegative values in',
     $      ' decreasing order, else 1/ulp',
     $      / ' 5 = | U - Upartial | / ( M ulp )',
     $      / ' 6 = | VT - VTpartial | / ( N ulp )',
     $      / ' 7 = | S - Spartial | / ( min(M,N) ulp |S| )',
     $      / ' ZGESDD: ', /
     $      ' 8 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ',
     $      / ' 9 = | I - U**T U | / ( M ulp ) ',
     $      / '10 = | I - VT VT**T | / ( N ulp ) ',
     $      / '11 = 0 if S contains min(M,N) nonnegative values in',
     $      ' decreasing order, else 1/ulp',
     $      / '12 = | U - Upartial | / ( M ulp )',
     $      / '13 = | VT - VTpartial | / ( N ulp )',
     $      / '14 = | S - Spartial | / ( min(M,N) ulp |S| )',
     $      / ' ZGESVJ: ', /
     $      / '15 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ',
     $      / '16 = | I - U**T U | / ( M ulp ) ',
     $      / '17 = | I - VT VT**T | / ( N ulp ) ',
     $      / '18 = 0 if S contains min(M,N) nonnegative values in',
     $      ' decreasing order, else 1/ulp',
     $      / ' ZGESJV: ', /
     $      / '19 = | A - U diag(S) VT | / ( |A| max(M,N) ulp )',
     $      / '20 = | I - U**T U | / ( M ulp ) ',
     $      / '21 = | I - VT VT**T | / ( N ulp ) ',
     $      / '22 = 0 if S contains min(M,N) nonnegative values in',
     $      ' decreasing order, else 1/ulp',
     $      / ' ZGESVDX(V,V,A): ', /
     $        '23 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ',
     $      / '24 = | I - U**T U | / ( M ulp ) ',
     $      / '25 = | I - VT VT**T | / ( N ulp ) ',
     $      / '26 = 0 if S contains min(M,N) nonnegative values in',
     $      ' decreasing order, else 1/ulp',
     $      / '27 = | U - Upartial | / ( M ulp )',
     $      / '28 = | VT - VTpartial | / ( N ulp )',
     $      / '29 = | S - Spartial | / ( min(M,N) ulp |S| )',
     $      / ' ZGESVDX(V,V,I): ',
     $      / '30 = | U**T A VT**T - diag(S) | / ( |A| max(M,N) ulp )',
     $      / '31 = | I - U**T U | / ( M ulp ) ',
     $      / '32 = | I - VT VT**T | / ( N ulp ) ',
     $      / ' ZGESVDX(V,V,V) ',
     $      / '33 = | U**T A VT**T - diag(S) | / ( |A| max(M,N) ulp )',
     $      / '34 = | I - U**T U | / ( M ulp ) ',
     $      / '35 = | I - VT VT**T | / ( N ulp ) ',
     $      ' ZGESVDQ(H,N,N,A,A',
     $      / '36 = | A - U diag(S) VT | / ( |A| max(M,N) ulp ) ',
     $      / '37 = | I - U**T U | / ( M ulp ) ',
     $      / '38 = | I - VT VT**T | / ( N ulp ) ',
     $      / '39 = 0 if S contains min(M,N) nonnegative values in',
     $      ' decreasing order, else 1/ulp',
     $      / / )
 9997 FORMAT( ' M=', I5, ', N=', I5, ', type ', I1, ', IWS=', I1,
     $      ', seed=', 4( I4, ',' ), ' test(', I2, ')=', G11.4 )
 9996 FORMAT( ' ZDRVBD: ', A, ' returned INFO=', I6, '.', / 9X, 'M=',
     $      I6, ', N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ),
     $      I5, ')' )
 9995 FORMAT( ' ZDRVBD: ', A, ' returned INFO=', I6, '.', / 9X, 'M=',
     $      I6, ', N=', I6, ', JTYPE=', I6, ', LSWORK=', I6, / 9X,
     $      'ISEED=(', 3( I5, ',' ), I5, ')' )

      RETURN

      // End of ZDRVBD

      }
