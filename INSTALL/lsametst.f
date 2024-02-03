      int                I1, I2;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ICHAR
      // ..
      // .. Executable Statements ..


      // Determine the character set.

      I1 = ICHAR( 'A' )
      I2 = ICHAR( 'a' )
      if ( I2-I1 == 32 ) {
         WRITE( *, * ) ' ASCII String    set';
      } else {
         WRITE( *, * ) ' Non-ASCII String    set, IOFF should be ',I2-I1;
      }

      // Test LSAME.

      IF( .NOT. LSAME( 'A', 'A' ) ) WRITE( *, 9999 )'A', 'A'       IF( .NOT. LSAME( 'A', 'a' ) ) WRITE( *, 9999 )'A', 'a'       IF( .NOT. LSAME( 'a', 'A' ) ) WRITE( *, 9999 )'a', 'A'       IF( .NOT. LSAME( 'a', 'a' ) ) WRITE( *, 9999 )'a', 'a'       IF( LSAME( 'A', 'B' ) ) WRITE( *, 9998 )'A', 'B'       IF( LSAME( 'A', 'b' ) ) WRITE( *, 9998 )'A', 'b'       IF( LSAME( 'a', 'B' ) ) WRITE( *, 9998 )'a', 'B'       IF( LSAME( 'a', 'b' ) ) WRITE( *, 9998 )'a', 'b'       IF( LSAME( 'O', '/' ) ) WRITE( *, 9998 )'O', '/'       IF( LSAME( '/', 'O' ) ) WRITE( *, 9998 )'/', 'O'       IF( LSAME( 'o', '/' ) ) WRITE( *, 9998 )'o', '/'       IF( LSAME( '/', 'o' ) ) WRITE( *, 9998 )'/', 'o'
      WRITE( *, * )' Tests completed'

 9999 FORMAT( ' *** Error:  LSAME( ', A1, ', ', A1, ') is false ' )
 9998 FORMAT( ' *** Error:  LSAME( ', A1, ', ', A1, ') is true ' )
      }
