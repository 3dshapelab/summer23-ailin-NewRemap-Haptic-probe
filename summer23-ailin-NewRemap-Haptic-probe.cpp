// this script aims to use probe adjustment task to measure the gain/strength of either texture or disparity as a depth cue
#include "stdafx.h"

void initStimulus(double stimulusDepth) {

	if (testing_texture_vs_disparity) {
		// testing texture, build vertices for curved surface
		buildVertices_textureMap(stimulusDepth);
	}
	else {
		/*
		// testing disparity, build random dots
		dot_container.clear();
		//dot_container = buildRandomDots(stimulusDepth);
		dot_container = buildRandomDots_Lattice(stimulusDepth, dot_num_per_col, dot_num_per_col, dot_jitter_max_scale);
		*/	
		buildVertices_randomDot(stimulusDepth);
		//buildVertices_panels_wavy(stimulusDepth);
		buildVertices_panels_jigsaw(stimulusDepth);
		
	}

	display_distance_jittered = display_distance + jitter_z;

}


void buildVertices_panels_wavy(double shapeDepth){

	vertContainer_Lcontour.clear();
	vertContainer_Rcontour.clear();

	int nr_seg = 1, nr_wave_perSeg = 6, nr_unit_perWave = 8, nr_step_unit = 10; // number of segments need to be an even number
	int del_unit_prev = 0, del_unit_next = 0, nr_unit_current = nr_unit_perWave;

	double step_y = stimulus_height / (nr_seg * nr_wave_perSeg * nr_unit_perWave * nr_step_unit);

	int rand_portions = 20; //how many portions we want - will creat (1 + rand_portions) different values
	double rand_min = .5, rand_max = 2;
	double seg_x_off_L, seg_x_off_R;

	double x_L = -stimulus_width / 2, x_R = stimulus_width / 2; // 
	double y = -stimulus_height / 2, y0 = -stimulus_height / 2;
	double tempy_height = ((double)(1.0/6.0)) * (stimulus_height / nr_seg);
	double dummy_sign = -1;

	double z = 0;
	vertContainer_Lcontour.push_back(Vector3d(x_L, y, z));
	vertContainer_Rcontour.push_back(Vector3d(x_R, y, z));

	// Left
	for(int i_wave = 0; i_wave < nr_wave_perSeg; i_wave++){

		if(i_wave < (nr_wave_perSeg - 1)){
			del_unit_next = rand() % 5 - 2;
		}else{
			del_unit_next = 0;
		}

		nr_unit_current = nr_unit_perWave - del_unit_prev + del_unit_next;
		tempy_height = nr_unit_current * nr_step_unit * step_y;


		if(i_wave % 2 == 0){

			dummy_sign = -1 * dummy_sign;	
			seg_x_off_L = rand_min + (rand_max - rand_min) * ( rand() % (rand_portions + 1) )/((double)(rand_portions));
			//seg_x_off_R = rand_min + (rand_max - rand_min) * ( rand() % (rand_portions + 1) )/((double)(rand_portions));


			for (int i_step = 0; i_step < nr_unit_current * nr_step_unit; i_step++){

				y = y + step_y;	
				x_L = -stimulus_width / 2 - dummy_sign * seg_x_off_L * sin(0.5 * M_PI * (y - y0) / tempy_height);
				//x_R = stimulus_width / 2 + dummy_sign * seg_x_off_R * sin(0.5 * M_PI * (y - y0) / tempy_height);	
				z = shapeDepth * cos(M_PI * y / stimulus_height);

				vertContainer_Lcontour.push_back(Vector3d(x_L, y, z));
				//vertContainer_Rcontour.push_back(Vector3d(x_R, y, z));
			
			}
			
		}else{

			for (int i_step = 0; i_step < nr_unit_current * nr_step_unit; i_step++){

				y = y + step_y;	
				x_L = -stimulus_width / 2 - dummy_sign * seg_x_off_L * cos(0.5 * M_PI * (y - y0) / tempy_height);
				//x_R = stimulus_width / 2 + dummy_sign * seg_x_off_R * cos(0.5 * M_PI * (y - y0) / tempy_height);	
				z = shapeDepth * cos(M_PI * y / stimulus_height);

				vertContainer_Lcontour.push_back(Vector3d(x_L, y, z));
				//vertContainer_Rcontour.push_back(Vector3d(x_R, y, z));
			
			}
			
		}

		y0 = y;
		del_unit_prev = del_unit_next;

	}

	// Right

	del_unit_prev = 0;
	y = -stimulus_height / 2, y0 = -stimulus_height / 2;
	dummy_sign = -1;
	for(int i_wave = 0; i_wave < nr_wave_perSeg; i_wave++){

		if(i_wave < (nr_wave_perSeg - 1)){
			del_unit_next = rand() % 5 - 2;
		}else{
			del_unit_next = 0;
		}

		nr_unit_current = nr_unit_perWave - del_unit_prev + del_unit_next;
		tempy_height = nr_unit_current * nr_step_unit * step_y;

		if(i_wave % 2 == 0){

			dummy_sign = -1 * dummy_sign;	
			//seg_x_off_L = rand_min + (rand_max - rand_min) * ( rand() % (rand_portions + 1) )/((double)(rand_portions));
			seg_x_off_R = rand_min + (rand_max - rand_min) * ( rand() % (rand_portions + 1) )/((double)(rand_portions));

			for (int i_step = 0; i_step < nr_unit_current * nr_step_unit; i_step++){

				y = y + step_y;	
				//x_L = -stimulus_width / 2 - dummy_sign * seg_x_off_L * sin(0.5 * M_PI * (y - y0) / tempy_height);
				x_R = stimulus_width / 2 + dummy_sign * seg_x_off_R * sin(0.5 * M_PI * (y - y0) / tempy_height);	
				z = shapeDepth * cos(M_PI * y / stimulus_height);

				//vertContainer_Lcontour.push_back(Vector3d(x_L, y, z));
				vertContainer_Rcontour.push_back(Vector3d(x_R, y, z));
			
			}
			
		}else{

			for (int i_step = 0; i_step < nr_unit_current * nr_step_unit; i_step++){

				y = y + step_y;	
				//x_L = -stimulus_width / 2 - dummy_sign * seg_x_off_L * cos(0.5 * M_PI * (y - y0) / tempy_height);
				x_R = stimulus_width / 2 + dummy_sign * seg_x_off_R * cos(0.5 * M_PI * (y - y0) / tempy_height);	
				z = shapeDepth * cos(M_PI * y / stimulus_height);

				//vertContainer_Lcontour.push_back(Vector3d(x_L, y, z));
				vertContainer_Rcontour.push_back(Vector3d(x_R, y, z));
			
			}
			
		}

		y0 = y;
		del_unit_prev = del_unit_next;

	}

}

void buildVertices_panels_jigsaw(double shapeDepth){

	vertContainer_Rcontour.clear();
	vertContainer_Lcontour.clear();
	//double width_separation = dispDist/(dispDist - shapeDepth) * stimulus_height * ratio_width_height;

	int nr_seg = 4, nr_sub = 10, nr_sub_firsthalf = 0; // number of segments need to be an even number
	double step_y = stimulus_height / (nr_seg * nr_sub);

	int rand_portions = 20; //how many steps we want - we will have (1 + rand_portions) different values
	double rand_min = .5, rand_max = 4;
	double seg_x_off_head = -(rand_min + (rand_max - rand_min) * ( rand() % (rand_portions + 1) )/((double)(rand_portions)) );
	double seg_x_off_end, seg_x_off_half, step_x_off;

	double x = stimulus_width / 2 + seg_x_off_head; // 
	double y = -stimulus_height / 2;
	double z = 0;
	vertContainer_Rcontour.push_back(Vector3d(x, y, z));


	for(int i_seg = 0; i_seg < nr_seg; i_seg++){

		seg_x_off_half = rand_min + (rand_max - rand_min) * ( rand() % (rand_portions + 1) )/((double)(rand_portions));
		seg_x_off_end = -(rand_min + (rand_max - rand_min) * ( rand() % (rand_portions + 1) )/((double)(rand_portions)) );
				
		nr_sub_firsthalf = nr_sub / 2 + (rand() % 5 - 2);

		step_x_off = (seg_x_off_half - seg_x_off_head) / nr_sub_firsthalf;

		for (int i_sub = 0; i_sub < nr_sub_firsthalf; i_sub++){

			x = x + step_x_off;
			y = y + step_y;			
			z = shapeDepth * cos(M_PI * y / stimulus_height);

			vertContainer_Rcontour.push_back(Vector3d(x, y, z)); 
			//vertContainer_Lcontour.push_back(Vector3d(-stimulus_width / 2, y, z));
		
		}

		step_x_off = (seg_x_off_end - seg_x_off_half) / (nr_sub - nr_sub_firsthalf);

		for (int i_sub = 0; i_sub < (nr_sub - nr_sub_firsthalf); i_sub++){

			x = x + step_x_off;
			y = y + step_y;			
			z = shapeDepth * cos(M_PI * y / stimulus_height);

			vertContainer_Rcontour.push_back(Vector3d(x, y, z));
			//vertContainer_Lcontour.push_back(Vector3d(-stimulus_width / 2, y, z));
		
		}

		seg_x_off_head = seg_x_off_end;
		
	}


	seg_x_off_head = (rand_min + (rand_max - rand_min) * ( rand() % (rand_portions + 1) )/((double)(rand_portions)) );

	x = -stimulus_width / 2 + seg_x_off_head; // 
	y = -stimulus_height / 2;
	z = 0;

	vertContainer_Lcontour.push_back(Vector3d(x, y, z));

	for(int i_seg = 0; i_seg < nr_seg; i_seg++){

		seg_x_off_end = rand_min + (rand_max - rand_min) * ( rand() % (rand_portions + 1) )/((double)(rand_portions));
		seg_x_off_half = -(rand_min + (rand_max - rand_min) * ( rand() % (rand_portions + 1) )/((double)(rand_portions)) );
				
		nr_sub_firsthalf = nr_sub / 2 + (rand() % 5 - 2);

		step_x_off = (seg_x_off_half - seg_x_off_head) / nr_sub_firsthalf;

		for (int i_sub = 0; i_sub < nr_sub_firsthalf; i_sub++){

			x = x + step_x_off;
			y = y + step_y;			
			z = shapeDepth * cos(M_PI * y / stimulus_height);

			vertContainer_Lcontour.push_back(Vector3d(x, y, z)); 
		
		}

		step_x_off = (seg_x_off_end - seg_x_off_half) / (nr_sub - nr_sub_firsthalf);

		for (int i_sub = 0; i_sub < (nr_sub - nr_sub_firsthalf); i_sub++){

			x = x + step_x_off;
			y = y + step_y;			
			z = shapeDepth * cos(M_PI * y / stimulus_height);

			vertContainer_Lcontour.push_back(Vector3d(x, y, z));
		
		}

		seg_x_off_head = seg_x_off_end;
		
	}
}



void drawVertices_panels(){

	int n = int(vertContainer_Rcontour.size());;
	float panel_width = 40;
	float panel_height_extra = 20;

	if( n > 0 ){

		glTranslated(0, 0, 2);

		glBegin(GL_QUAD_STRIP);

		glVertex3f(vertContainer_Lcontour.at(0)[0],		vertContainer_Lcontour.at(0)[1] - panel_height_extra,			vertContainer_Lcontour.at(0)[2]); //0
		glVertex3f(vertContainer_Lcontour.at(0)[0] -  panel_width,					vertContainer_Lcontour.at(0)[1] - panel_height_extra,			vertContainer_Lcontour.at(0)[2]); //1

		for (int i = 0; i < n; i++)
		{	
			glVertex3f(vertContainer_Lcontour.at(i)[0],		vertContainer_Lcontour.at(i)[1],			vertContainer_Lcontour.at(i)[2]); //0
			glVertex3f(vertContainer_Lcontour.at(i)[0] - panel_width,					vertContainer_Lcontour.at(i)[1],			vertContainer_Lcontour.at(i)[2]); //1

		}	

		glVertex3f(vertContainer_Lcontour.at(n-1)[0],		vertContainer_Lcontour.at(n-1)[1] + panel_height_extra,			vertContainer_Lcontour.at(n-1)[2]); //0
		glVertex3f(vertContainer_Lcontour.at(n-1)[0] - panel_width,					vertContainer_Lcontour.at(n-1)[1] + panel_height_extra,			vertContainer_Lcontour.at(n-1)[2]); //1

		glEnd();

		// Right panel
		glBegin(GL_QUAD_STRIP);

		glVertex3f(vertContainer_Rcontour.at(0)[0] + panel_width,		vertContainer_Rcontour.at(0)[1] - panel_height_extra,			vertContainer_Rcontour.at(0)[2]); //0
		glVertex3f(vertContainer_Rcontour.at(0)[0],					vertContainer_Rcontour.at(0)[1] - panel_height_extra,			vertContainer_Rcontour.at(0)[2]); //1

		for (int i = 0; i < n; i++)
		{	
			glVertex3f(vertContainer_Rcontour.at(i)[0] + panel_width,		vertContainer_Rcontour.at(i)[1],			vertContainer_Rcontour.at(i)[2]); //0
			glVertex3f(vertContainer_Rcontour.at(i)[0],					vertContainer_Rcontour.at(i)[1],			vertContainer_Rcontour.at(i)[2]); //1

		}	

		glVertex3f(vertContainer_Rcontour.at(n-1)[0] + panel_width,		vertContainer_Rcontour.at(n-1)[1] + panel_height_extra,			vertContainer_Rcontour.at(n-1)[2]); //0
		glVertex3f(vertContainer_Rcontour.at(n-1)[0],					vertContainer_Rcontour.at(n-1)[1] + panel_height_extra,			vertContainer_Rcontour.at(n-1)[2]); //1

		glEnd();
	}

}

void buildVertices_randomDot(double shapeDepth) {
	
	int nr_pts_height = 101;
	int nr_pts_width = (nr_pts_height - 1) * ratio_width_height_long + 1;

	double long_width = stimulus_height * ratio_width_height_long;

	double step_size_width = long_width / (double)(nr_pts_width - 1);
	double step_size_height = stimulus_height / (double)(nr_pts_height - 1);

	GLuint i_ind = 0;

	double x, y, z;
	//double normal_x, normal_y, normal_z;

	vertices_vec_new.clear();
	colors_vec_new.clear();
	indices_draw_triangle_vec_new.clear();


	for (int j = 0; j < nr_pts_height; j++) {  // 

		y = -stimulus_height / 2 + j * step_size_height;
		z = shapeDepth * cos(M_PI * y / stimulus_height);
		
		//normal_x = 0;
		//normal_y = shapeDepth * sin(M_PI * y / stimulus_height) * M_PI / stimulus_height;
		//normal_z = 1;
		

		for (int i = 0; i < nr_pts_width; i++) { //

			x = -long_width / 2 + i * step_size_width;
	
			// using vector
			vertices_vec_new.push_back(x);
			vertices_vec_new.push_back(y);
			vertices_vec_new.push_back(z);

			//normals_vec.push_back(normal_x);
			//normals_vec.push_back(normal_y);
			//normals_vec.push_back(normal_z);

			if(rand() % 8 == 1){
				colors_vec_new.push_back(0);
				colors_vec_new.push_back(0);
				colors_vec_new.push_back(0);
			}else{
				colors_vec_new.push_back(0.9);
				colors_vec_new.push_back(0);
				colors_vec_new.push_back(0);
			}


			// construct the triangle indices to be drawn
			if (i < nr_pts_width - 1 && j < nr_pts_height - 1) {

				indices_draw_triangle_vec_new.push_back(i_ind);
				indices_draw_triangle_vec_new.push_back(i_ind + 1);
				indices_draw_triangle_vec_new.push_back(i_ind + nr_pts_width);

				indices_draw_triangle_vec_new.push_back(i_ind + nr_pts_width);
				indices_draw_triangle_vec_new.push_back(i_ind + 1);
				indices_draw_triangle_vec_new.push_back(i_ind + nr_pts_width + 1);
				//ind = ind + 6;
			}

			i_ind++;

		}

	}
	
}

void drawVertices_randomDot(double dispDist, double shapeDepth){

	glPushMatrix();

	glLoadIdentity();
	glTranslated(0, 0, dispDist - shapeDepth);

	// enable matrices for use in drawing below
	glEnable(GL_POLYGON_SMOOTH);
	glEnable(GL_BLEND);
	//glEnable(GL_NORMALIZE); //so we don't need to normalize our normal for surfaces

	// activate and specify pointer to vertex array
	glEnableClientState(GL_VERTEX_ARRAY);
	//glEnableClientState(GL_TEXTURE_COORD_ARRAY);
	//glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);

	glVertexPointer(3, GL_FLOAT, 0, &vertices_vec_new[0]);
	//glTexCoordPointer(2, GL_FLOAT, 0, &texcoors_vec[0]);
	//glNormalPointer(GL_FLOAT, 0, &normals_vec[0]); //
	glColorPointer(3, GL_FLOAT, 0, &colors_vec_new[0]);
	glDrawElements(GL_TRIANGLES, indices_draw_triangle_vec_new.size(), GL_UNSIGNED_INT, &indices_draw_triangle_vec_new[0]);

	// deactivate vertex arrays after drawing
	glDisableClientState(GL_VERTEX_ARRAY);
	//glDisableClientState(GL_TEXTURE_COORD_ARRAY);
	//glDisableClientState(GL_NORMAL_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);


	switch(panel_state){
	case no_aperture:
		break;

	case black_aperture:
		glColor3f(0.0f, 0.0f, 0.0f);		
		drawVertices_panels();
		//such that the separation has the same visual angle of the stimulusheight
		break;

	case red_aperture:
		glColor3f(0.5f, 0.0f, 0.0f);		
		drawVertices_panels();
		break;


	}

	glPopMatrix();
}

void buildVertices_textureMap(double shapeDepth) {
	
	double step_size_width = stimulus_width / (double)(nr_points - 1);
	double step_size_height = stimulus_height / (double)(nr_points - 1);

	GLuint i_ind = 0;

	double x, y, z, y_prev, z_prev, u, v;
	y_prev = -stimulus_height / 2; 
	z_prev = 0; 
	double total_distance_y = 0; //tracks the distance along y/z axis, approximate the "diameter" of the ellipse
	double normal_x, normal_y, normal_z;

	vertices_vec.clear();
	colors_vec.clear();
	texcoors_vec.clear();
	normals_vec.clear();
	indices_draw_triangle_vec.clear();


	for (int j = 0; j < nr_points; j++) {  // 

		y = -stimulus_height / 2 + j * step_size_height;
		z = shapeDepth * cos(M_PI * y / stimulus_height);


		total_distance_y = total_distance_y + sqrt(pow(y - y_prev, 2) + pow(z - z_prev, 2));
		v = total_distance_y / normalizer_to_uv + v_offset; //v coordinate
		
		normal_x = 0;
		normal_y = shapeDepth * sin(M_PI * y / stimulus_height) * M_PI / stimulus_height;
		normal_z = 1;
		

		for (int i = 0; i < nr_points; i++) { //

			x = -stimulus_width / 2 + i * step_size_width;
			u = (x + stimulus_width / 2) / normalizer_to_uv + u_offset;; //u coordinate. 

		/* 
	
		step 1: build the meshgrid using vertices, 


			(8)---(9)---(10)--(11)
			 |     |     |     |
			 |     |     |     |
			 |     |     |     |
			 |     |     |     |
			(4)---(5)---(6)---(7)
			 |     |     |     |
			 |     |     |     |
			 |     |     |     |
			 |     |     |     |
			(0)---(1)---(2)---(3)

		
			going over each vertex: store the (x,y,z) to vertices array or vector, (u,v) to texcoors array or vector, (1,0,0) to colors array or vector,
			if light source or shading is involved, normals are needed

*/			
	
			// using vector
			vertices_vec.push_back(x);
			vertices_vec.push_back(y);
			vertices_vec.push_back(z);

			colors_vec.push_back(1);
			colors_vec.push_back(0);
			colors_vec.push_back(0);

			texcoors_vec.push_back(u);
			texcoors_vec.push_back(v);

			normals_vec.push_back(normal_x);
			normals_vec.push_back(normal_y);
			normals_vec.push_back(normal_z);
	/*
		step 2: create an array/vector that store how the triangles should be drawn
		
		The array/vector is one dimensional, but is groupped by unit of 3, meaning every three elements form a triangle

		for example, if indices_draw_triangle is like [0 1 4 4 1 5 1 2 5 5 2 6 ...], 
		then it draws triangles with indices {0 1 4}, {4 1 5}, {1 2 5}, {5 2 6}...

		(8)---(9)---(10)--(11)
		 |\    |\    |\    |
		 | \   | \   | \   |
		 |  \  |  \  |  \  |
		 |   \ |   \ |   \ |
		(4)---(5)---(6)---(7)
		 |\    |\    |\    |  
		 | \   | \   | \   | 
		 |  \  |  \  |  \  |
		 |   \ |   \ |   \ |
		(0)---(1)---(2)---(3)



	 triangle 1: 0 1 4;   triangle 2: 4 1 5 
	 triangle 3: 1 2 5;   triangle 2: 5 2 6 ...
*/

			// construct the triangle indices to be drawn
			if (i < nr_points - 1 && j < nr_points - 1) {

				indices_draw_triangle_vec.push_back(i_ind);
				indices_draw_triangle_vec.push_back(i_ind + 1);
				indices_draw_triangle_vec.push_back(i_ind + nr_points);

				indices_draw_triangle_vec.push_back(i_ind + nr_points);
				indices_draw_triangle_vec.push_back(i_ind + 1);
				indices_draw_triangle_vec.push_back(i_ind + nr_points + 1);
				//ind = ind + 6;
			}

			i_ind++;

		}
				
		y_prev = y; z_prev = z;
	}
	
}

void drawVertices_textureMap(int texNum, double dispDist, double shapeDepth) {

	glShadeModel(GL_SMOOTH); // enable Smooth Shading
	glEnable(GL_LIGHTING); // enable lighting
	glEnable(GL_LIGHT1); 
	glEnable(GL_NORMALIZE); //so we don't need to normalize our normal for surfaces

	// Light source parameters
	GLfloat LightAmbient[]= { amb_intensity, 0.0f, 0.0f, 1.0f }; // non-directional & overall light (r,g,b,alpha): dark part
	GLfloat LightDiffuse[]= { 1 - amb_intensity, 0.0f, 0.0f, 1.0f }; // light created by the light source (directional light; r,g,b,alpha): bright part
	GLfloat LightPosition[]= { 0.0f, 1.f, 0.8f, 0.0f }; // Light Position (x, y, z, 1.0f); if w==0, directional; if w==1, positional lights. Attenuation can be applied only to the positional light 

	//glPushMatrix();
	//glLoadIdentity();

	//setting the light
	glLightfv(GL_LIGHT1, GL_AMBIENT, LightAmbient); //setup the ambient light
	glLightfv(GL_LIGHT1, GL_DIFFUSE, LightDiffuse); //setup the diffuse light
	glLightfv(GL_LIGHT1, GL_POSITION,LightPosition); //position the light

	glPushMatrix();

	glLoadIdentity();
	glTranslated(0, 0, dispDist - shapeDepth);

	// enable matrices for use in drawing below
	//glEnable(GL_LIGHTING);
	glEnable(GL_POLYGON_SMOOTH);
	glEnable(GL_BLEND);
	glEnable(GL_TEXTURE_2D);
	//glEnable(GL_NORMALIZE); //so we don't need to normalize our normal for surfaces

	// bind the texture

	glBindTexture(GL_TEXTURE_2D, loaded_textures[texNum]);

	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	// activate and specify pointer to vertex array
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_TEXTURE_COORD_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);

	glVertexPointer(3, GL_FLOAT, 0, &vertices_vec[0]);
	glTexCoordPointer(2, GL_FLOAT, 0, &texcoors_vec[0]);
	glNormalPointer(GL_FLOAT, 0, &normals_vec[0]); //
	glColorPointer(3, GL_FLOAT, 0, &colors_vec[0]);
	glDrawElements(GL_TRIANGLES, indices_draw_triangle_vec.size(), GL_UNSIGNED_INT, &indices_draw_triangle_vec[0]);

	glPopMatrix();

	// deactivate vertex arrays after drawing
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_TEXTURE_COORD_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);

	glDisable(GL_LIGHTING); 
}



std::vector<Vector3d> buildRandomDots(double shapeDepth)
{
	std::vector<Vector3d> dotContainer;

	for (int dots_placed = 0; dots_placed < dot_number; dots_placed++)
	{

		double x = (rand() % 101)/100.0 * stimulus_width - stimulus_width / 2;
		double y = (rand() % 101)/100.0 * stimulus_height - stimulus_height / 2;

		double z = shapeDepth * cos(M_PI * y / stimulus_height);


		dotContainer.push_back(Vector3d(x, y, z));
	}

	return dotContainer;
}


std::vector<Vector3d> buildRandomDots_Lattice(double shapeDepth, int dotNumPerRow, int dotNumPerCol, double dotJitterMax_Scale)
{
	std::vector<Vector3d> dotContainer;

	int dot_counter_top = 0, dot_counter_bottom = 0;
	int N_count = 2;

	double step_size_x = stimulus_width / (double)(dotNumPerRow - 1);
	double step_size_y = stimulus_height / (double)(dotNumPerCol - 1);

	double dot_jitter_x_max = step_size_x * dotJitterMax_Scale;
	double dot_jitter_y_max = step_size_y * dotJitterMax_Scale;

	double x,y,z;

	for(int i_y = 0; i_y < dotNumPerCol; i_y++){
			for(int i_x = 0; i_x < dotNumPerRow; i_x++){
				x = i_x * step_size_x - stimulus_width / 2 +
					(rand() % 17) / 16.0 * dot_jitter_x_max - dot_jitter_x_max / 2;

				y = i_y * step_size_y - stimulus_height / 2 +
					(rand() % 17) / 16.0 * dot_jitter_y_max - dot_jitter_y_max / 2;

				if(y < - stimulus_height / 2){
					y = -stimulus_height /2;
					z = 0;
					dot_counter_bottom++;

					if(dot_counter_bottom % N_count == 0){
						dotContainer.push_back(Vector3d(x, y, z));
					}

				}else if(y > stimulus_height / 2){
					y = stimulus_height / 2;
					z = 0;
					dot_counter_top++;

					if(dot_counter_top % N_count == 0){
						dotContainer.push_back(Vector3d(x, y, z));
					}
				}else{

					z = shapeDepth * cos(M_PI * y / stimulus_height);
					dotContainer.push_back(Vector3d(x, y, z));
				}
				
			}
	}

	dot_number = (int)(dotContainer.size()/ratio_width_height);


	return dotContainer;
}


void drawRandomDots(std::vector<Vector3d> dotContainer, double dispDist, double shapeDepth) {
	
	glPushMatrix();

	glLoadIdentity();
	glTranslated(0, 0, dispDist - shapeDepth);

	glColor3f(dot_color_R, 0.0f, 0.0f);

	for (int i = 0; i < int(dotContainer.size()); i++)
	{
		Vector3d dot_vector = dotContainer.at(i);
		double x_axis = dot_vector[0];
		double y_axis = dot_vector[1];
		double z_axis = dot_vector[2];

		glPushMatrix();
		glTranslated(x_axis, y_axis, z_axis);

		double dot_size = ratio_dotSize_distance * abs(display_distance + z_axis);
		glutSolidSphere(dot_size, dot_sides, dot_sides);
		glPopMatrix();
	}

	glPopMatrix();

	switch(panel_state){
		case no_aperture:
			break;

		case black_aperture:
			glPushMatrix();
			glLoadIdentity();
			glTranslated(0, 0, dispDist + 2);
			glColor3f(0.0f, 0.0f, 0.0f);
			
			drawBlockingPanels((dispDist/(dispDist - shapeDepth)) * stimulus_height); 
			//such that the separation has the same visual angle of the stimulusheight
			glPopMatrix();
			break;

		case red_aperture:
			glPushMatrix();
			glLoadIdentity();
			glTranslated(0, 0, dispDist + 2);
			glColor3f(0.4f, 0.0f, 0.0f);
			
			drawBlockingPanels((dispDist/(dispDist - shapeDepth)) * stimulus_height);
			glPopMatrix();
			break;


	}

	
}


void drawBlockingPanels(double pandelSeparation){

	double panel_w = 20, panel_h = 80;

	// left panel
	glDisable(GL_TEXTURE_2D);

	glBegin(GL_QUADS);
	glVertex3f(-pandelSeparation / 2 - panel_w,   panel_h / 2,  0.0f);
	glVertex3f(-pandelSeparation / 2,             panel_h / 2,  0.0f);
	glVertex3f(-pandelSeparation / 2,            -panel_h / 2,  0.0f);
	glVertex3f(-pandelSeparation / 2 - panel_w,  -panel_h / 2,  0.0f);
	glEnd();

	// right panel
	glBegin(GL_QUADS);
	glVertex3f(pandelSeparation / 2,              panel_h / 2,  0.0f);
	glVertex3f(pandelSeparation / 2 + panel_w,    panel_h / 2,  0.0f);
	glVertex3f(pandelSeparation / 2 + panel_w,   -panel_h / 2,  0.0f);
	glVertex3f(pandelSeparation / 2,             -panel_h / 2,  0.0f);
	glEnd();
}

void drawFixation(double dispDist) {
	// draws a small fixation cross at the center of the display
	glDisable(GL_TEXTURE_2D);
	glColor3f(1.0f, 0.0f, 0.0f);
	glLineWidth(2.f);

	glPushMatrix();
	glLoadIdentity();
	glTranslated(0, 0, dispDist);
	double cross_length = 5;
	glBegin(GL_LINES);
	glVertex3d(cross_length / 2, 0, 0);
	glVertex3d(-cross_length / 2, 0, 0);
	glVertex3d(0, -cross_length / 2. , 0);
	glVertex3d(0, cross_length / 2. , 0);
	glEnd();

	glPopMatrix();

}

// This function seems to be used to shut down the system after use
void shutdown(){
	cout << "shutting down" << endl;
	responseFile.close(); // close this object
	if(resetScreen_betweenRuns){
		homeEverything(5000,4500);}
	cleanup();
	exit(0);
}
void cleanup() 
{
// Stop the optotrak
    optotrak->stopCollection();
    delete optotrak;
}
void initProjectionScreen(double _focalDist, const Affine3d &_transformation, bool synchronous)
{
	focalDistance = _focalDist;	
    screen.setWidthHeight(SCREEN_WIDE_SIZE, SCREEN_WIDE_SIZE*SCREEN_HEIGHT/SCREEN_WIDTH);
    screen.setOffset(alignmentX,alignmentY);
    screen.setFocalDistance(_focalDist);
    screen.transform(_transformation);
    cam.init(screen);
	if ( synchronous )
		moveScreenAbsolute(_focalDist,homeFocalDistance,4500);
	else
		moveScreenAbsoluteAsynchronous(_focalDist,homeFocalDistance,4500);
}
// Initialize Optotrak for use in the experiment
void initOptotrak()
{
    optotrak=new Optotrak2(); //intiailize the Optotrak object
    optotrak->setTranslation(calibration);

	//define Optotrak-specific variables
    int numMarkers=22;
    float frameRate=85.0f;
    float markerFreq=4600.0f;
    float dutyCycle=0.4f;
    float voltage = 7.0f;

	// run the intiailization method for the Optotrak, checking to see if ever (if == 0) and catch the error if so
    if ( optotrak->init("C:/cncsvisiondata/camerafiles/Aligned20111014",numMarkers, frameRate, markerFreq, dutyCycle,voltage) != 0)
    {   cerr << "Something during Optotrak initialization failed, press ENTER to continue. A error log has been generated, look \"opto.err\" in this folder" << endl;
        cin.ignore(1E6,'\n');
        exit(0);
    }

    // Read 10 frames of coordinates and fill the markers vector
    for (int i=0; i<10; i++)
    {
        updateTheMarkers();
    }
}
// run a method to define a vector that holds marker positions 
void updateTheMarkers()
{
	optotrak->updateMarkers();
	markers = optotrak->getAllMarkers();

}
// Initialize motors for moving screen around
void initMotors()
{
	//specify the speed for (objects,screen)
		if(resetScreen_betweenRuns){
		homeEverything(5000,4500);}

}

// Method that initializes the openGL parameters needed for creating the stimuli. 
// seems like this is not changed for each experiment (maybe for different experimental setup eg monitor)
void initRendering()
{   

	glClearColor(0.0,0.0,0.0,1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	/* Set depth buffer clear value */
	glClearDepth(1.0);
	/* Enable depth test */
	glEnable(GL_DEPTH_TEST);
	/* Set depth function */
	glDepthFunc(GL_LEQUAL);
	// scommenta solo se vuoi attivare lo shading degli stimoli

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glLineWidth(1.5);

	/*
	if(testing_texture_vs_disparity){

		//texture-only
		glEnable(GL_MULTISAMPLE);
		glHint(GL_MULTISAMPLE_FILTER_HINT_NV, GL_NICEST);
		
	}else{
		// disparity-only, random dots
		glEnable(GL_LINE_SMOOTH);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
		glEnable(GL_POINT_SMOOTH);
		glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	}
*/	
		glEnable(GL_MULTISAMPLE);
		glHint(GL_MULTISAMPLE_FILTER_HINT_NV, GL_NICEST);
}


void initVariables()
{
	interoculardistance = atof(parameters.find("IOD").c_str());
	display_distance = str2num<double>(parameters.find("dispDepth"));
	
	switch(sessionNum){
		case 1:
			testing_texture_vs_disparity = false;
			left_eye_on = true;
			right_eye_on = true;
			break;

		case 2: 
			testing_texture_vs_disparity = true;
			left_eye_on = true;
			right_eye_on = false;
			break;

		case 3:
			testing_texture_vs_disparity = false;
			left_eye_on = true;
			right_eye_on = true;
			break;

		case 4: 
			testing_texture_vs_disparity = true;
			left_eye_on = false;
			right_eye_on = true;
			break;

	}

	monocular_display = testing_texture_vs_disparity;

	
	// eye coordinates
	if(testing_texture_vs_disparity){
		// if testing texture, then both eys are set to be eyeMiddle
		eyeRight = Vector3d(0, 0, 0);
		eyeLeft = Vector3d(0, 0, 0);
	}else{
		eyeRight = Vector3d(interoculardistance / 2, 0, 0);
		eyeLeft = Vector3d(-interoculardistance / 2, 0, 0);
	}

	
	stimulus_height = tan((DEG2RAD * visual_angle)/2) * 2 * (abs(display_distance));
	stimulus_width =  ratio_width_height * stimulus_height;

}

void initBlock()
{
	// initialize the trial matrix
	trial.init(parameters);
	trial.next();

	trialNum = 1;
}

// Initialize the streams, open the file and write to it
void initStreams()
{
	// Initialize the parameter file starting from the file parameters.txt, if the file does not exist, it tells you
	ifstream parametersFile;

	//parametersFileName = experiment_directory + "/parameters_summer22-ailin-hapticRemap-SingleDepthCue-probe.txt";

	parametersFile.open(parametersFileName.c_str());
	parameters.loadParameterFile(parametersFile);


	// Subject name
	subjectName = parameters.find("SubjectName");
	string session = parameters.find("Session");
	sessionNum = str2num<int>(parameters.find("Session"));

	if ((sessionNum < 1 ) || (sessionNum > 4 ))
	{
		string error_on_file_io = string("invalid session number");
		cerr << error_on_file_io << endl;
		MessageBox(NULL, (LPCSTR)"INVALID SESSION",NULL, NULL);
		shutdown();
	}

	string responseFileName = experiment_directory +"/"+ subjectName + "_s" + session + ".txt";
	// Principal streams files
	if (util::fileExists(experiment_directory +"/"+ subjectName + "_s" + session + ".txt") && subjectName != "junk")
	{
		string error_on_file_io = string("file already exists");
		cerr << error_on_file_io << endl;
		MessageBox(NULL, (LPCSTR)"FILE ALREADY EXISTS\n Please check the parameters file.",NULL, NULL);
		shutdown();
	}

	responseFile.open(responseFileName.c_str());
	responseFile << fixed << responseFile_headers << endl;
	
	//globalTimer.start();
}

void drawGLScene()
{
		glDrawBuffer(GL_BACK_RIGHT);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glClearColor(0.0,0.0,0.0,1.0);
		cam.setEye(eyeLeft);

		if(left_eye_on){
			drawInfo();
			drawStimulus();
		}

		
		// Draw right eye view
		glDrawBuffer(GL_BACK_LEFT);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glClearColor(0.0,0.0,0.0,1.0);
		cam.setEye(eyeRight); 

		if(right_eye_on){
			drawInfo();
			drawStimulus();
		}

	glutSwapBuffers();
}

void update(int value)
{
    glutPostRedisplay();
    glutTimerFunc(TIMER_MS, update, 0);
}

void drawProgressBar() {
	
	glPushMatrix();
	glLoadIdentity();
	glTranslated(0, 0, display_distance);

	glColor3f(0.2, 0.2, 0.2);
	glBegin(GL_LINE_LOOP);
	glVertex3f(-50, 5, 0);
	glVertex3f(50, 5, 0);
	glVertex3f(50, -5, 0);
	glVertex3f(-50, -5, 0);
	glEnd();

	glColor3f(0.1, 0.3, 0.1);
	glBegin(GL_POLYGON);
	glVertex3f(-50, 5, 0);
	glVertex3f(-50 + percentComplete, 5, 0);
	glVertex3f(-50 + percentComplete, -5, 0);
	glVertex3f(-50, -5, 0);
	glEnd();

	glPopMatrix();
}

void drawTimeBar(){

    glPushMatrix();
    glLoadIdentity();
    glTranslated(0, -150, display_distance);
	//glColor3fv(glRed);
	if(ElapsedTime < responseTime){
		glColor3f(0, 0.5, 0);
	}else{
		glColor3f(0.6, 0, 0);
		}
    glLineWidth(2.0);
    //slash 1
    glBegin(GL_LINE_LOOP);
    glVertex3d(-50,0,0);
    glVertex3d(50,0,0);
    glEnd();
    //slash 2

    glBegin(GL_LINE_LOOP);
    glVertex3d(-50 + ElapsedTime/(responseTime/100),5,0);
    glVertex3d(-50 + ElapsedTime/(responseTime/100),-5,0);
    glEnd();
    glPopMatrix();
    glLineWidth(1.0);
}


void drawInfo()
{
	// displays relevant information to the screen
	if ( visibleInfo )
	{
		glDisable(GL_COLOR_MATERIAL);
		glDisable(GL_BLEND);
		GLText text;
		if (gameMode)
			text.init(SCREEN_WIDTH, SCREEN_HEIGHT, glWhite, GLUT_BITMAP_HELVETICA_18);
		else
			text.init(640, 480, glWhite, GLUT_BITMAP_HELVETICA_12);

	
		text.enterTextInputMode();

		switch(current_stage){
			case stimulus_preview:
				glColor3fv(glRed);
				text.draw("Welcome! press + to start training");
				text.draw("# Name: " + subjectName);
				text.draw("# IOD: " + stringify<double>(interoculardistance));
				text.draw("# depth: " + stringify<double>(depth_test));
				text.draw("# left eye on: " + stringify<bool>(left_eye_on));
				text.draw("# rt eye on: " + stringify<bool>(right_eye_on));
				text.draw("# IOD: " + stringify<float>(dot_color_R));
				
				
					// check if mirror is calibrated
				
				text.draw("# !!!!Mirror Alignment = " + stringify<double>(mirrorAlignment));

				break;

			case trial_fixate:
			case trial_view:
			case trial_respond:

				glColor3fv(glRed);
				text.draw("# Name: " + subjectName);
				text.draw("# IOD: " + stringify<double>(interoculardistance));
				text.draw("# depth: " + stringify<double>(depth_test));
				text.draw("# trial: " + stringify<int>(trialNum));
				text.draw("# time: " + stringify<double>(ElapsedTime));
				text.draw("# current stage: " + stringify<int>(current_stage));

									// check if mirror is calibrated
				if (abs(mirrorAlignment - 45.0) > 0.2)
					text.draw("# !!!!Mirror Alignment = " + stringify<double>(mirrorAlignment));
				break;

			case break_time:
				if(training && trialNum >= trainNum_cap)
					text.draw("Pls call the experimenter!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
				else
					text.draw("Break time! Press + to continue");
				break;

			case exp_completed:
				glColor3fv(glWhite);
				text.draw("The experiment is over. Thank you! :)");
				break;
		}
		text.leaveTextInputMode();
		glEnable(GL_COLOR_MATERIAL);
		glEnable(GL_BLEND);
	}
}


// Funzione che gestisce il ridimensionamento della finestra
void handleResize(int w, int h)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glViewport(0,0,SCREEN_WIDTH, SCREEN_HEIGHT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

}

void drawStimulus()
{
	//enum Stages { stimulus_preview, prep_trial, trial_view, trial_respond, trial_distractor, wait_next_trial, break_time, exp_completed };
	switch(current_stage){

		case stimulus_preview:

			if (testing_texture_vs_disparity) {
				// testing texture, build vertices for curved surface
				drawVertices_textureMap(texnum, display_distance_jittered, depth_test);
				
			}
			else {
				// testing disparity, build random dots
				//drawRandomDots(dot_container, display_distance_jittered, depth_test);
				drawVertices_randomDot(display_distance_jittered, depth_test);
			}

			break;

		case trial_fixate:
			drawFixation(display_distance_jittered);
			break;

		case trial_view:
			if (testing_texture_vs_disparity) {
				// testing texture, build vertices for curved surface
				drawVertices_textureMap(texnum, display_distance_jittered, depth_test);
			}
			else {
				// testing disparity, build random dots
				drawRandomDots(dot_container, display_distance_jittered, depth_test);
			}
			break;

		case trial_respond:

			if (testing_texture_vs_disparity) {
				// testing texture, build vertices for curved surface
				drawVertices_textureMap(texnum, display_distance_jittered, depth_test);
			}
			else {
				// testing disparity, build random dots
				drawRandomDots(dot_container, display_distance_jittered, depth_test);
			}

			drawProbe(display_distance_jittered + probe_jitter_z, probe_depth);
			break;

		case break_time:
			drawProgressBar();
			break;



	}
}
/*
double getZ(double theHeight, double theDepth, double currentY){
	
	double currentZ;

	switch(current_shape){

		case Ridge:		

		if(theDepth < theHeight / 2 )
		{
			double R = (pow(theDepth, 2) + pow((theHeight / 2.), 2)) / (2 * theDepth);

			if (abs(currentY) >= theHeight / 2) {
				currentZ = 0;
			}
			else {
				currentZ = sqrt(pow(R, 2) - pow(currentY, 2)) - R + theDepth;
			}
		}
		else{

			if (abs(currentY) >= theHeight / 2) {
				currentZ = 0;
			}
			else {
				currentZ = theDepth * sqrt(1 - pow(currentY / (theHeight / 2.), 2));
			}
		}

		break;

		case Gaussian:
			//double sig_height_ratio = 0.14;
			currentZ = theDepth * exp(-pow(currentY,2)/(2 * pow(gauss_sig_height_ratio * theHeight, 2)));

			break;

		case Cosine:
			//double phase_y = M_PI * currentY / (theHeight / 2.);
			currentZ = (theDepth / 2) * cos(M_PI * currentY / (theHeight / 2.)) + (theDepth / 2);

			break;

		case CosineRidge:
			currentZ = theDepth * cos(M_PI * currentY / theHeight) ;
			break;
	}

	return currentZ;

}
*/


void updateProbe(double probeHeight, double probeDepth){

	double step_size = probeHeight / (double)(nr_points - 1);

	double y; 

	for (int j = 0; j < nr_points; j++) { 

		y = -probeHeight / 2 + j * step_size;
		probe_y_array[j] = y;
		probe_x_array[j] = probeDepth * cos(M_PI * y / probeHeight) ;
	}

}

void drawProbe(double dispDist, double probeDepth){

	glDisable(GL_TEXTURE_2D);
	//glColor3f(1.0f, 0.0f, 0.0f);

	glPushMatrix();
	glLoadIdentity();
	glTranslated(0, probe_location_offset, dispDist);

	double prev_x, prev_y, x, y;
	x = probe_x_array[0];
	y = probe_y_array[0];

	glBegin(GL_LINES);
	glLineWidth(1.f);
	glColor3f(0.7f, 0.0f, 0.0f);
	for(int i = 1; i < nr_points; i++){

		prev_x = x;
		prev_y = y;

		x = probe_x_array[i];
		y = probe_y_array[i];


		glVertex3d(prev_x, prev_y, 0);
		glVertex3d(x, y, 0);
	}

	glEnd();

	glBegin(GL_LINES);
	glLineWidth(2.f);
	glColor3f(0.8f, 0.0f, 0.0f);
	glVertex3d(0, -40, 0);
	glVertex3d(0,  40, 0);
	glVertex3d(probeDepth, -40, 0);
	glVertex3d(probeDepth,  40, 0);
	glEnd();


	glPopMatrix();
}

void initTrial()
{
	current_stage = prep_trial;
	initProjectionScreen(display_distance);
	trial_timer.reset();
	trial_timer.start();
	if(training){
		depth_test = rand() % depth_training_range + depth_training_min ;

		//shapeID = rand() % 2 + 2;
		//current_shape = (shapeTypes)(shapeID);

	}
	else {
		depth_test = trial.getCurrent()["testDepths"];
		//normalizer_to_uv = trial.getCurrent()["textNormalizer"];

		//shapeID =  trial.getCurrent()["testShape"];
		//current_shape = (shapeTypes)(shapeID);

	}
	if(monocular_display){
		jitter_z = (rand() % 101 - 50)/10.0;
		probe_jitter_z = - depth_test;
	}else{
		jitter_z = (rand() % 101 - 50)/10.0;
		probe_jitter_z = -1.0 * (rand() % (int)(depth_test));
	}
	// build the stimulus and choose texture
	initStimulus(depth_test);
	texnum = rand() % 50 + 1;
	normalizer_to_uv = normalizer_to_uv_base + ((rand() % 21) - 10); // from -15 to 15
	
	//probe init
	probe_depth_init = 1. + rand() % (int) (1.2 * depth_training_range);
	probe_depth = probe_depth_init;
	updateProbe(stimulus_height, probe_depth);
	
		
	current_stage = trial_fixate;

	trial_timer.reset();				
	trial_timer.start();
	ElapsedTime = 0;

	// reset ElapsedTime variables

	//trial_timer.start(); // start the timer when key is pressed

	
	
}

void onlineTrial(){

	switch(current_stage){


	case trial_fixate:

		if (ElapsedTime > fixateTime) {
			current_stage = trial_view;
		}
		break;

	case trial_view:

		if (ElapsedTime > fixateTime + viewTime) {
			current_stage = trial_respond;
		}
		break;
	}

}

void advanceTrial()
{
	//subjName\tIOD\tblockN\ttrialN\tdisplayDistance\tvisualAngle\tclyHorizontal\ttexnum\ttextNomralizer\ttestDepth\tprobeStart\tprobeDepth\ttime
	if(training){
		if(trialNum < trainNum_cap){
			beepOk(4);
			trialNum++;
			initTrial();
		}else{
			beepOk(2);
			trialNum++;
			current_stage = break_time;
			visibleInfo = true;			
		}

	}else{
		//subjName\tIOD\tblockN\ttrialN\tdisplayDistance\tvisualAngle\tshapeID\ttexnum\ttextNomralizer\ttestDepth\tprobeDepthInit\tprobeDepth\tRT
		responseFile << fixed <<
		subjectName << "\t" <<		
		interoculardistance << "\t" <<
		blkNum << "\t" <<
		trialNum << "\t" <<
		display_distance << "\t" <<
		visual_angle << "\t" <<
		stimulus_height << "\t" <<
		(int)current_shape << "\t" <<
		testing_texture_vs_disparity << "\t" <<
		!testing_texture_vs_disparity<< "\t" <<
		left_eye_on << "\t" <<
		right_eye_on << "\t" <<
		texnum << "\t" <<
		normalizer_to_uv << "\t" <<
		depth_test << "\t" <<
		probe_depth_init << "\t" <<
		probe_depth << "\t" << 
		dot_number << "\t" << 
		visAngle_dot << "\t" << 
		dot_color_R << "\t" << 
		dot_jitter_max_scale << "\t" << 
		ElapsedTime  << endl;	

		if (!trial.isEmpty()){

			if(trialNum % 100 == 0){
				beepOk(4);				
				percentComplete = trialNum /(totalTrNum/100.0);
				trial.next();				
				current_stage = break_time;
				visibleInfo = true;

			}else{
				beepOk(4);
				trialNum++;
				trial.next();
				initTrial();
							}
			
		}else{
			beepOk(1);
			responseFile.close();
			visibleInfo=true;
			current_stage = exp_completed;
		}

	}

}

void setLightingPar(){

	glShadeModel(GL_SMOOTH); // enable Smooth Shading

	// Light source parameters
	GLfloat LightAmbient[]= { amb_intensity, 0.0f, 0.0f, 1.0f }; // non-directional & overall light (r,g,b,alpha): dark part
	GLfloat LightDiffuse[]= { 1 - amb_intensity, 0.0f, 0.0f, 1.0f }; // light created by the light source (directional light; r,g,b,alpha): bright part
	GLfloat LightPosition[]= { 0.0f, 1.f, 0.8f, 0.0f }; // Light Position (x, y, z, 1.0f); if w==0, directional; if w==1, positional lights. Attenuation can be applied only to the positional light 

	glPushMatrix();
	glLoadIdentity();

	//setting the light
	glLightfv(GL_LIGHT1, GL_AMBIENT, LightAmbient); //setup the ambient light
	glLightfv(GL_LIGHT1, GL_DIFFUSE, LightDiffuse); //setup the diffuse light
	glLightfv(GL_LIGHT1, GL_POSITION,LightPosition); //position the light

	glPopMatrix();

}


void handleKeypress(unsigned char key, int x, int y)
{   
	switch (key){ // key presses that work regardless of the stage
	
		case 27:	//corrisponde al tasto ESC
			shutdown();		
		break;

		case 'i':
			visibleInfo=!visibleInfo;
		break;

		case 'o':
			panel_state = panelStates((panel_state + 1)%3);
			break;

		case '1':
			switch(current_stage){

				case stimulus_preview:

					if(depth_test > depth_inc)
						depth_test = depth_test - depth_inc;

					initStimulus(depth_test);
					break;

				case trial_respond:
						if(probe_depth > 1)
							probe_depth = probe_depth - 1;

						updateProbe(stimulus_height, probe_depth);
						break;

				case prep_trial:
				case trial_fixate:
				case trial_view:
				case break_time:
				case exp_completed:
					beepOk(3);
					break;

			}
			break;

		case '2':
			switch(current_stage){

				case stimulus_preview:

					depth_test = depth_test + depth_inc;
					initStimulus(depth_test);
					break;

				case trial_respond:
					
					probe_depth = probe_depth + 1;
					updateProbe(stimulus_height, probe_depth);
					break;

				case prep_trial:
				case trial_fixate:
				case trial_view:
				case break_time:
				case exp_completed:
					beepOk(3);
					break;

			}
			break;


		case '+':
			
			switch(current_stage){
				case stimulus_preview:
					
					beepOk(5);
					visibleInfo = false;
					initBlock();
					initTrial();
					
					//initStimulus(depth_test);
				break;

				case trial_respond:
					trial_timer.stop();
					advanceTrial();
				break;

				case break_time:
					beepOk(5);
					visibleInfo = false;
					trialNum++;
					initTrial();
				break;		
			}
			break;

		case 'T':
		case 't':
			if(current_stage!= stimulus_preview){
				if(training){
					training = false;
					beepOk(6);
					trialNum = 0;
					current_stage = break_time;
				}
			}
			break;
			/////////////////////////
			// adjusting disparity //
			/////////////////////////


		case '7':
			if(dot_color_R > .2)
				dot_color_R = dot_color_R - .05;
			initStimulus(depth_test);
			break;

		case '8':
			dot_color_R = dot_color_R + .05;
			initStimulus(depth_test);
			break;

		case '6':
			if(dot_num_per_col > 3)
				dot_num_per_col = dot_num_per_col - 1;
			initStimulus(depth_test);
			break;

		case '9':
			dot_num_per_col = dot_num_per_col + 1;
			initStimulus(depth_test);
			break;

		///////////////////////////////////////////////
		// key presses for adjusting stimulus preview

		case 'A':
		case 'a':
			if(current_stage == stimulus_preview){
				testing_texture_vs_disparity = !testing_texture_vs_disparity;
				initStimulus(depth_test);
			}
			break;

		case 'M':
		case 'm':
			if(current_stage == stimulus_preview){
				//jitter_z = jitter_z - 2;
				//display_distance_jittered = display_distance + jitter_z;
				current_shape = shapeTypes((current_shape + 1) % 4);
				initStimulus(depth_test);
			}
			break;


	}

}

/***** SOUNDS *****/
void beepOk(int tone)
{

	switch(tone)
	{

	case 1: //high pitch beep
		PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\beep-8_lowpass.wav", NULL, SND_FILENAME | SND_ASYNC);
	break;

	case 4: //mellow and good for trials
		PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\beep-440-pluck.wav", NULL, SND_FILENAME | SND_ASYNC);
	break;

	case 2: //reject like buzz
		PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\beep-10.wav", NULL, SND_FILENAME | SND_ASYNC);
	break;

	case 3: //reject short
		PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\beep-reject.wav", NULL, SND_FILENAME | SND_ASYNC);
	break;

	case 5: //mellow and good for trials
		PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\beep-highBubblePop.wav", NULL, SND_FILENAME | SND_ASYNC);
	break;

	case 6: //mellow and good for trials
		PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\beep-lowBubblePop.wav", NULL, SND_FILENAME | SND_ASYNC);
	break;

	case 15:
		PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\beep-rising.wav",
			NULL, SND_FILENAME | SND_ASYNC);
	break;

	case 16:
		PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\beep-falling.wav",
			NULL, SND_FILENAME | SND_ASYNC);
	break;

	case 17:
	PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\beep-440-pluck-5below.wav",
		NULL, SND_FILENAME | SND_ASYNC);
	break;

	case 18: // light click
	PlaySound((LPCSTR) "C:\\cygwin\\home\\visionlab\\workspace\\cncsvision\\data\\beep\\beep-click3MS.wav",
		NULL, SND_FILENAME | SND_ASYNC);
	break;

	}
	return;
}
void idle()
{
	onlineTrial();
	online_apparatus_alignment();
	ElapsedTime = trial_timer.getElapsedTimeInMilliSec();
}

/*** Online operations ***/
void online_apparatus_alignment()
{
	updateTheMarkers();
    // mirror alignment check
    mirrorAlignment = asin(
        abs((markers.at(mirror1).p.z() - markers.at(mirror2).p.z())) /
        sqrt(
        pow(markers.at(mirror1).p.x() - markers.at(mirror2).p.x(), 2) +
        pow(markers.at(mirror1).p.z() - markers.at(mirror2).p.z(), 2)
        )
        ) * 180 / M_PI;

    // screen Y alignment check
    screenAlignmentY = asin(
		abs((markers.at(screen1).p.y() - markers.at(screen3).p.y())) /
		sqrt(
		pow(markers.at(screen1).p.x() - markers.at(screen3).p.x(), 2) +
		pow(markers.at(screen1).p.y() - markers.at(screen3).p.y(), 2)
		)
		) * 180 / M_PI;

    // screen Z alignment check
    screenAlignmentZ = asin(
		abs(markers.at(screen1).p.z() - markers.at(screen2).p.z()) /
		sqrt(
		pow(markers.at(screen1).p.x() - markers.at(screen2).p.x(), 2) +
		pow(markers.at(screen1).p.z() - markers.at(screen2).p.z(), 2)
		)
		) * 180 / M_PI *
		abs(markers.at(screen1).p.x() - markers.at(screen2).p.x()) /
		(markers.at(screen1).p.x() - markers.at(screen2).p.x());

}

int LoadGLTextures()  // Load PNG And Convert To Textures
{

	for (int i = 1; i <= 50; i++) {
		std::stringstream ss;
		ss << i;

		string texturePath = "fall16-ailin-stimulusTest/0/polkadots" + ss.str() + ".png";
		 loaded_textures[i] = SOIL_load_OGL_texture
		(
			texturePath.c_str(),
			SOIL_LOAD_AUTO,
			SOIL_CREATE_NEW_ID,
			SOIL_FLAG_MULTIPLY_ALPHA
		);
	}


	return true; // Return Success
}

// this is run at compilation because it's titled 'main'
int main(int argc, char*argv[])  
{
	//functions from cncsvision packages
	mathcommon::randomizeStart();
		
	// initializing glut (to use OpenGL)
	glutInit(&argc, argv);
	/*
	if(testing_texture_vs_disparity){
		glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_STEREO | GLUT_MULTISAMPLE);
	}else{
		glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STEREO);
	}
*/

	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_STEREO | GLUT_MULTISAMPLE);

	glutGameModeString("1024x768:32@85"); //resolution  
	glutEnterGameMode();
	glutFullScreen();
	
	// initializes optotrak and velmex motors
	initOptotrak();
	initMotors();
	
	initRendering(); // initializes the openGL parameters needed for creating the stimuli
	
	initStreams(); // streams as in files for writing data

	initVariables();

	LoadGLTextures();

	initStimulus(depth_test);

	initProjectionScreen(display_distance);

	// glut callback, OpenGL functions that are infinite loops to constantly run 

	glutDisplayFunc(drawGLScene); // keep drawing the stimuli

	glutKeyboardFunc(handleKeypress); // check for keypress

	glutReshapeFunc(handleResize);

	glutIdleFunc(idle);

	glutTimerFunc(TIMER_MS, update, 0);

	glutSetCursor(GLUT_CURSOR_NONE);

	//boost::thread initVariablesThread(&initVariables); 

	glutMainLoop();

	//cleanup();
	return 0;
}
