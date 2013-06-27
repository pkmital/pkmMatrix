/*
 *  pkmImage.h
 *  pkmSourceSeparation_iPhone
 *
 *  Created by Mr. Magoo on 7/3/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#pragma once

#include <Accelerate/Accelerate.h>
#include "ofImage.h"
#include "pkmMatrix.h"
class pkmImage
{
public:
	pkmImage()
	{
		data_scaled = NULL;
	}
	~pkmImage()
	{
		if (bAllocated) {
			free(data_scaled);
		}
	}
	
	void allocate(int rows, int cols)
	{
		data_scaled = (float *)malloc(sizeof(float) * rows * cols);
		pixels = (char *)malloc(sizeof(char) * rows * cols);
		image.allocate(cols, rows, OF_IMAGE_GRAYSCALE);
		bAllocated = true;
		r = rows;
		c = cols;
		len = r*c;
	}
	
	void convertFloatToOFImage(float *data, bool normalize = true)
	{
		if (normalize) {
			float max, min;
			vDSP_maxv(data, 1, &max, len);
			vDSP_minv(data, 1, &min, len);
			
			min = -min;
			vDSP_vsadd(data, 1, &min, data_scaled, 1, len);
			min = -min;
			float scale = 255.0 / (max - min);
			vDSP_vsmul(data_scaled, 1, &scale, data_scaled, 1, len);
			vDSP_vfix8(data_scaled, 1, pixels, 1, len);
			//image.setFromPixels((unsigned char *)pixels, c, r, OF_IMAGE_GRAYSCALE);
			memcpy(image.getPixels(), pixels, len*sizeof(unsigned char));
			image.update();
		}
		else {
			float scale = 255.0;
			vDSP_vsmul(data, 1, &scale, data_scaled, 1, len);
			vDSP_vfix8(data_scaled, 1, pixels, 1, len);
			memcpy(image.getPixels(), pixels, len*sizeof(unsigned char));
			//image.setFromPixels((unsigned char *)pixels, c, r, OF_IMAGE_GRAYSCALE);
			image.update();
		}
	}
	
	int len;
	int r,c;
	char *pixels;
	float *data_scaled;
	ofImage image;
	bool bAllocated;
};