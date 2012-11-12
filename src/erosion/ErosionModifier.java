package erosion;

import java.util.Random;

import src.Array2D;
import src.IModifier;

public class ErosionModifier implements IModifier<Array2D>
{
	private static final float FLOAT_EPSILON = 1E-9f;
	
	public static final class Point
	{
		public float x, y;
		
		public Point(float x, float y)
		{
			this.x = x;
			this.y = y;
		}
	}
	
	public static final class PointMap
	{
		private final Point[] map;
		public final int width;
		
		public PointMap(int width, int height)
		{
			this.width = width;

			map = new Point[width * height];
			
			for (int i = 0; i < map.length; ++i)
				map[i] = new Point(0, 0);
		}
		
		public final Point get(int x, int y)
		{
			return map[y * width + x];
		}
	}
	
	private final void deposit(Array2D hmap, PointMap erosion, int xi, int zi, float xf, float zf, float ds)
	{
		depositAt(hmap, erosion, xi,     zi,     ds * (1f - xf) * (1f - zf));
		depositAt(hmap, erosion, xi + 1, zi,     ds * xf        * (1f - zf));
		depositAt(hmap, erosion, xi,     zi + 1, ds * (1f - xf) * zf);
		depositAt(hmap, erosion, xi + 1, zi + 1, ds * xf        * zf);
		/*
		depositAt(xi  , zi  , 0.25f) \
		depositAt(xi+1, zi  , 0.25f) \
		depositAt(xi  , zi+1, 0.25f) \
		depositAt(xi+1, zi+1, 0.25f) \
		 */
	}
	
	private final void depositAt(Array2D hmap, PointMap erosion, int x, int z, float delta)
	{
	    erosion.get(x, z).y += delta;
	    hmap.set(x, z, Math.min(1f, hmap.get(x, z) + delta));
	    //params.deposit(surface[HMAP_INDEX(X, Z)], delta);
	}
	
	private final void erode(Array2D hmap, PointMap erosion, int x, int z, float delta)
	{
		hmap.set(x, z, Math.max(0f, hmap.get(x, z) - delta));
		Point e = erosion.get(x, z);
		
		float r = e.x, d = e.y;
		
		if (delta <= d) 
		{
			d -= delta;
		}
		else
		{
			r += delta - d;
			d = 0;
		}
		
		e.x = r;
		e.y = d;
		//scolor=params.erode(surface[HMAP_INDEX(X, Z)], s, delta);
	}
	
	private final Random random;
	private final int iterations;
	private final float Kq, waterEvapRate, soilErodeRate, soilDepositRate, inertia, minSlope, gravity;
	private final IErosionStyle erosionStyle;
	
	public ErosionModifier(
			Random random,
			int iterations,
			float Kq,
			float waterEvapRate,
			float soilErodeRate,
			float soilDepositRate,
			float inertia,
			float minSlope,
			float gravity,
			IErosionStyle erosionStyle)
	{
		this.random = random;
		this.iterations = iterations;
		this.Kq = Kq;
		this.waterEvapRate = waterEvapRate;
		this.soilErodeRate = soilErodeRate;
		this.soilDepositRate = soilDepositRate;
		this.inertia = inertia;
		this.minSlope = minSlope;
		this.gravity = gravity;
		this.erosionStyle = erosionStyle;
	}
	
	public void modify(Array2D heightMap)
	{
		PointMap erosion = new PointMap(heightMap.width, heightMap.height);
		
		//float flt = 0;

		//long t0 = get_ref_time();

		int maxPathLength = (int)(Math.hypot(heightMap.width, heightMap.height) * 10d);
		long longPaths = 0, randomDirs = 0, sumLen = 0;
		
		
		erosionStyle.setup(heightMap);

		for (int iter = 0; iter < iterations; ++iter)
		{
			//if ((iter & 0x3FFF) == 0 && iter != 0)
				//show_splash("Calculating erosion", (iter + 0.5f) / iterations);
			
			IErosionStyle.Point start = erosionStyle.makeStartPoint();
			int xi = start.x;
			int zi = start.y;
			
			float xp = xi, zp = zi;
			float xf = 0, zf = 0;
			
			float h = heightMap.get(xi, zi);
			
			float carried_sediment = erosionStyle.makeInitialSediment();
			float velocity = erosionStyle.makeInitialVelocity();
			float water = 1f;
			//vec4f scolor = zero4f();
			
			float h00 = h;
			float h10 = heightMap.get(xi + 1, zi  );
			float h01 = heightMap.get(xi,     zi + 1);
			float h11 = heightMap.get(xi + 1, zi + 1);
			
			float dx = 0, dz = 0;
				
			int numMoves = 0;
	    	  
			for (; numMoves < maxPathLength; ++numMoves)
			{
				if (xi <= 0 || xi >= heightMap.width  - 2
				 || zi <= 0 || zi >= heightMap.height - 2)
				{
					break;
				}
				
			    // calc gradient
			    float gx = h00 + h01 - h10 - h11;
			    float gz = h00 + h10 - h01 - h11;
			    //== better interpolated gradient?
		
			    // calc next pos
			    dx = (dx - gx) * inertia + gx;
			    dz = (dz - gz) * inertia + gz;

			    float dl = (float)Math.hypot(dx, dz);
			    if (dl <= FLOAT_EPSILON)
			    {
			        // pick random dir
			        float a = random.nextFloat() * (float)Math.PI * 2f;
			        dx = (float)Math.cos(a);
			        dz = (float)Math.sin(a);
			        ++randomDirs;
			    }
			    else
			    {
			    	dx /= dl;
			    	dz /= dl;
			    }

				float xp_next = xp + dx;
				float zp_next = zp + dz;
				
				// sample next height
				int xi_next = (int)xp_next;
				int zi_next = (int)zp_next;
				float xf_next = xp_next - xi_next;
				float zf_next = zp_next - zi_next;
				
				float nh00 = heightMap.get(xi_next,     zi_next  );
				float nh10 = heightMap.get(xi_next + 1, zi_next  );
				float nh01 = heightMap.get(xi_next,     zi_next + 1);
				float nh11 = heightMap.get(xi_next + 1, zi_next + 1);
				
				float h_next = (nh00 * (1 - xf_next) + nh10 * xf_next) * (1 - zf_next)
						     + (nh01 * (1 - xf_next) + nh11 * xf_next) * zf_next;


				// if higher than current, try to deposit sediment up to neighbor height
			      
			    if (h_next >= h)
			    {
			    	float delta_sediment = (h_next - h) + 0.00001f;

			    	if (delta_sediment >= carried_sediment)
					{
					    // deposit all sediment and stop
						delta_sediment = carried_sediment;
						deposit(heightMap, erosion, xi, zi, xf, zf, carried_sediment);
						h += delta_sediment;
						carried_sediment = 0;
						break;
					}
		
			        deposit(heightMap, erosion, xi, zi, xf, zf, delta_sediment);
			        h += delta_sediment;
			        carried_sediment -= delta_sediment;
			        velocity = 0;
			    }

	      		// compute transport capacity
	      		float dh = h - h_next;
	      		float slope = dh;
	      		//float slope=dh / sqrtf(dh * dh + 1);

	      		float next_sediment = Math.max(slope, minSlope) * velocity * water * Kq;

	      		// deposit/erode (don't erode more than dh)
		        float delta_sediment = carried_sediment - next_sediment;
		        if (delta_sediment >= 0)
		        {
	    	  		// deposit
			        delta_sediment *= soilDepositRate;
			        //ds = Math.min(ds, 1.0f);
		
			        deposit(heightMap, erosion, xi, zi, xf, zf, delta_sediment);
			        dh += delta_sediment;
			        carried_sediment -= delta_sediment;
		        }
		        else
		        {
			        // erode
			        delta_sediment *= -soilErodeRate;
			        delta_sediment = Math.min(delta_sediment, dh * 0.99f);

					for (int z = zi - 1; z <= zi + 2; ++z)
					{
						float zo = z - zp;
						float zo2 = zo * zo;
							
						for (int x = xi - 1; x <= xi + 2; ++x)
						{
							float xo = x - xp;
							
							float W = 1f - (xo * xo + zo2) * 0.25f;
							
							if (W <= 0)
								continue;
							
							W *= 0.1591549430918953f;
						
							erode(heightMap, erosion, x, z, delta_sediment * W);
						}
					}
		        /*
		            ERODE(xi  , zi  , (1-xf)*(1-zf))
		            ERODE(xi+1, zi  ,    xf *(1-zf))
		            ERODE(xi  , zi+1, (1-xf)*   zf )
		            ERODE(xi+1, zi+1,    xf *   zf )*/
		
			        dh -= delta_sediment;
		
			        carried_sediment += delta_sediment;
	      		}

			  	// move to the neighbor
				velocity = (float)Math.sqrt(velocity * velocity + gravity * dh);
				water *= 1f - waterEvapRate;
				
				xp = xp_next; zp = zp_next;
				xi = xi_next; zi = zi_next;
				xf = xf_next; zf = zf_next;
				
				h = h_next;
				h00 = nh00;
				h10 = nh10;
				h01 = nh01;
				h11 = nh11;
			}

			if (numMoves >= maxPathLength)
			{
				System.out.println("droplet "+iter+" path is too long!");
				++longPaths;
			}
			
			sumLen += numMoves;
		}
		
		/*long t1 = get_ref_time();
		System.out.println(String.format("computed %7d erosion droplets in %6u ms, %.0f droplets/s",
		iterations, get_time_msec(t1 - t0), (double)iterations / get_time_sec(t1 - t0)));
		
		System.out.println(String.format("  %.2f average path length, %I64u long paths cut, %I64u random directions picked",
		(double)sumLen / iterations, longPaths, randomDirs));*/
	}
}
