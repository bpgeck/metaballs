
# include "my_viewer.h"

# include <sigogl/ws_run.h>

int main ( int argc, char** argv )
{
	MyViewer* v = new MyViewer(-1, -1, 1280, 720, "My SIG Application");
	v->cmd ( WsViewer::VCmdAxis );

	v->view_all ();
	v->show ();

	ws_run ();
	return 1;
}
