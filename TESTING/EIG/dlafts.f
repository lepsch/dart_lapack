      SUBROUTINE DLAFTS( TYPE, M, N, IMAT, NTESTS, RESULT, ISEED, THRESH, IOUNIT, IE )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             TYPE;
      int                IE, IMAT, IOUNIT, M, N, NTESTS
      double             THRESH;
*     ..
*     .. Array Arguments ..
      int                ISEED( 4 )
      double             RESULT( * );
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      int                K
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAHD2
*     ..
*     .. Executable Statements ..
*
      IF( M.EQ.N ) THEN
*
*     Output for square matrices:
*
         DO 10 K = 1, NTESTS
            IF( RESULT( K ).GE.THRESH ) THEN
*
*           If this is the first test to fail, call DLAHD2
*           to print a header to the data file.
*
               IF( IE.EQ.0 ) CALL DLAHD2( IOUNIT, TYPE )
               IE = IE + 1
               IF( RESULT( K ).LT.10000.0D0 ) THEN
                  WRITE( IOUNIT, FMT = 9999 )N, IMAT, ISEED, K, RESULT( K )
 9999             FORMAT( ' Matrix order=', I5, ', type=', I2,
     $                  ', seed=', 4( I4, ',' ), ' result ', I3, ' is',
     $                  0P, F8.2 )
               ELSE
                  WRITE( IOUNIT, FMT = 9998 )N, IMAT, ISEED, K, RESULT( K )
 9998             FORMAT( ' Matrix order=', I5, ', type=', I2,
     $                  ', seed=', 4( I4, ',' ), ' result ', I3, ' is',
     $                  1P, D10.3 )
               END IF
            END IF
   10    CONTINUE
      ELSE
*
*     Output for rectangular matrices
*
         DO 20 K = 1, NTESTS
            IF( RESULT( K ).GE.THRESH ) THEN
*
*              If this is the first test to fail, call DLAHD2
*              to print a header to the data file.
*
               IF( IE.EQ.0 ) CALL DLAHD2( IOUNIT, TYPE )
               IE = IE + 1
               IF( RESULT( K ).LT.10000.0D0 ) THEN
                  WRITE( IOUNIT, FMT = 9997 )M, N, IMAT, ISEED, K, RESULT( K )
 9997             FORMAT( 1X, I5, ' x', I5, ' matrix, type=', I2, ', s',
     $                  'eed=', 3( I4, ',' ), I4, ': result ', I3,
     $                  ' is', 0P, F8.2 )
               ELSE
                  WRITE( IOUNIT, FMT = 9996 )M, N, IMAT, ISEED, K, RESULT( K )
 9996             FORMAT( 1X, I5, ' x', I5, ' matrix, type=', I2, ', s',
     $                  'eed=', 3( I4, ',' ), I4, ': result ', I3,
     $                  ' is', 1P, D10.3 )
               END IF
            END IF
   20    CONTINUE
*
      END IF
      RETURN
*
*     End of DLAFTS
*
      END
