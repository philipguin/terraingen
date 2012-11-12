package erosion;

import src.Array2D;
import src.IDerivator;


public class WindElevationDerivator implements IDerivator<Array2D>
{
	private final float initialElevation, windDescentRate, minWindHeight;
	
	public WindElevationDerivator(float initialElevation, float windDescentRate, float minWindHeight)
	{
		this.initialElevation = initialElevation;
		this.windDescentRate = windDescentRate;
		this.minWindHeight = minWindHeight;
	}
	
	public Array2D derive(Array2D heightMap)
	{
		Array2D windMap = new Array2D(heightMap.width, heightMap.height);
		
		double angle = 40f;//random.nextDouble() * Math.PI * 2d;
		float dx = (float)Math.cos(angle);
		float dz = (float)Math.sin(angle);

		//x and z should be no greater than 1
		
		if (dx > dz)
		{
			dz /= Math.abs(dx);
			dx /= Math.abs(dx);
		}
		else
		{
			dx /= Math.abs(dz);
			dz /= Math.abs(dz);
		}
		
		float windDescentRate = this.windDescentRate * (float)Math.hypot(dx, dz);
		
		int h_align, v_align, v_start, v_end;
		
		if (dz >= 0f)
		{
			h_align = 0;
			v_start = 1;
			v_end = windMap.height;
		}
		else
		{
			h_align = windMap.height - 1;
			v_start = 0;
			v_end = windMap.height - 1;
		}

		if (dx >= 0f)
		{
			v_align = 0;
		}
		else
		{
			v_align = windMap.width - 1;
		}

		for (int i = 0; i < windMap.width;  ++i)
			iteration_v(windMap, heightMap, i, h_align, dx, dz, dz >= 0f ? 1 : -1, windDescentRate);
		
		for (int j = v_start; j < v_end; ++j)
			iteration_h(windMap, heightMap, v_align, j, dx, dz, dx >= 0f ? 1 : -1, windDescentRate);
		
		return windMap;
	}
	
	private final void iteration_h(Array2D windMap, Array2D heightMap, int x, int z, float dx, float dz, int xIt, float windDescentRate)
	{
		//Assume either dx or dz is 1f and the other is less.
		float z_precise = z;
		float elevation = initialElevation;
		
		int lastZ;

		dz /= Math.abs(dx);
		windDescentRate /= Math.abs(dx);
		
		iteration:
		while (true)
		{
			lastZ = z;

			while (z == lastZ)
			{
				elevation = Math.max(minWindHeight, elevation - windDescentRate);
				
				float height = heightMap.get(x, z);
				
				if (elevation < height)
					elevation = height;
				
				windMap.set(x, z, elevation - heightMap.get(x, z));

				x += xIt;
				z_precise += dz;
				z = (int)z_precise;
				
				if (x < 0 || x >= windMap.width || z < 0 || z >= windMap.height)
					break iteration;
			}
		}
	}
	
	private final void iteration_v(Array2D windMap, Array2D heightMap, int x, int z, float dx, float dz, int zIt, float windDescentRate)
	{
		//Assume either dx or dz is 1f and the other is less.
		float x_precise = x;
		float elevation = initialElevation;
		
		int lastX;
		
		dx /= Math.abs(dz);
		windDescentRate /= Math.abs(dz);
		
		iteration:
		while (true)
		{
			lastX = x;

			while (x == lastX)
			{
				elevation = Math.max(minWindHeight, elevation - windDescentRate);
				
				float height = heightMap.get(x, z);
				
				if (elevation < height)
					elevation = height;
				
				windMap.set(x, z, elevation - heightMap.get(x, z));

				z += zIt;
				x_precise += dx;
				x = (int)x_precise;

				if (x < 0 || x >= windMap.width || z < 0 || z >= windMap.height)
					break iteration;
			}
		}
	}
}
