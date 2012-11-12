package src;

import java.util.Random;

public class PerlinNoiseArray2DPopulator implements IPopulator<Array2D>
{
	private final Random random;
	private final int gran;
	private final float offset, range;
	
	public PerlinNoiseArray2DPopulator(Random random, int granularity, float min, float max)
	{
		this.random = random;
		this.gran = granularity;
		this.range = (max - min) * .5f;
		this.offset = min + range;
	}
	
	private final static float sqrtOf2 = 1.41421356F;//sqrtf(2.0);
	
	private static final float s_curve(float t)
	{
	    return t * t * (3f - 2f * t);
	}
	
	private int[] makeP()
	{       
	    int[] result = new int[256];
	    
	    for (int i = 0; i < 256; ++i)
	    	result[i] = random.nextInt() & 255;
	    
	    return result;
	}
	
	/*private static void initialize1D()
	{
	    g1D = new float[256];
	    
	    for (int i = 0; i < 256; ++i)
	        g1D[i] = rand_f() * 2f - 1f;
	}
	
	private static final float gradients1D(final int x)
	{
	    return g1D[x & 255]; 
	}
	
	public static float[] PerlinNoise1D(int width, int gran)
	{
	    initializeP();
	    initialize1D();
	    
	    //generate noise
	    float[] result = new float[width];
	    
	    float x, o,p;
	    int gI;
	    
	    for (int i = 0; i < width; ++i)
	    {
	        x = (float)(i%gran)/gran;
	        gI = i/gran;
	        
	        o = gradients1D(gI)   * x;
	        p = gradients1D(gI+1) * (x-1f);
	        
	        result[i] = o + s_curve(x)*(p-o);
	    }
	    
	    
	    return result;
	}*/
	
	private float[][] makeGradients()
	{
	    float[][] result = new float[256][];
	    float[] gradient;
	    double angle;
	    
	    for (int i = 0; i < 256; ++i)
	    {
	        gradient = result[i] = new float[2];
	        
	        angle = random.nextDouble() * 2d * Math.PI;
	        gradient[0] = (float)Math.sin(angle) * sqrtOf2 * range;
	        gradient[1] = (float)Math.cos(angle) * sqrtOf2 * range;
	    }
	    
	    return result;
	}
	
	private final float[] getGradient(float[][] gradients, int[] P, int x, int y)
	{
	    return gradients[(x + P[y & 255]) & 255]; 
	}
	
	public void populate(Array2D result)
	{
		int[] P = makeP();
	    float[][] gradients = makeGradients();
	    
	    int i, j;
	    
	    //generate noise
	    
	    float[] gradient;
	    float x,y, xInv,yInv, sx,sy, o,p,q,r;
	    int gI,gJ, gI1,gJ1;
	    long n = 0;
	    
	    int width = result.width;
	    int height = result.height;
	    
	    for (i = 0; i < width; ++i)
	    {
	        x = (float)(i % gran) / gran;
	        xInv = x - 1f;
	        gI = i / gran;
	        gI1 = gI + 1;
	        sx = s_curve(x);
	        
	        for (j = 0; j < height; ++j)
	        {
	            y = (float)(j % gran) / gran;
	            yInv = y - 1f;
	            gJ = j / gran;
	            gJ1 = gJ + 1;
	            sy = s_curve(y);
	            
	            gradient = getGradient(gradients, P, gI,  gJ);
	            o = gradient[0] * x    + gradient[1] * y;
	            
	            gradient = getGradient(gradients, P, gI1, gJ);
	            p = gradient[0] * xInv + gradient[1] * y;
	            
	            gradient = getGradient(gradients, P, gI,  gJ1);
	            q = gradient[0] * x    + gradient[1] * yInv;
	            
	            gradient = getGradient(gradients, P, gI1, gJ1);
	            r = gradient[0] * xInv + gradient[1] * yInv;
	                
	                
	            o =  o + sx * (p - o);
	            p =  q + sx * (r - q);
	            
	            result.set1D((int)n++, offset + o + sy * (p - o));
	        }
	    }
	}
	
	/*private static void initialize3D()
	{   
	    g3D = new float[256][];
	    float[] gradient;
	    float theta;
	    
	    for (int i = 0; i < 256; ++i)
	    {
	        gradient = g3D[i] = new float[3];
	        
	        theta = rand_f()*2f*3.14159265f;
	        gradient[0] = (float)Math.sin(theta);
	        gradient[1] = (float)Math.cos(theta);
	        
	        theta = rand_f()*2f*3.14159265f;
	        gradient[2] = (float)Math.sin(theta)*sqrtOf2;
	        
	        theta = (float)Math.cos(theta)*sqrtOf2;
	        gradient[0] *= theta;
	        gradient[1] *= theta;
	    }
	}
	
	private static final float[] gradients3D(final int x, final int y, final int z)
	{
	    return g3D[(x + p[(y + p[z & 255]) & 255]) & 255]; 
	}
	
	public static float[] PerlinNoise3D(int width, int height, int depth, int gran)
	{
	    initializeP();
	    initialize3D();
	    
	    int i, j, k;
	    
	    //generate noise
	    float[] result = new float[width*height*depth];
	    
	    float[] gradient;
	    float x,y,z, xInv,yInv,zInv, sx,sy,sz, o,p,q,r,s,t,u,v;
	    int gI,gJ,gK, gI1,gJ1,gK1;
	    long n = 0;
	    
	    for (i = 0; i < width; ++i)
	    {
	        x = (float)(i%gran)/gran;
	        xInv = x-1f;
	        gI = i/gran;
	        gI1 = gI+1;
	        sx = s_curve(x);
	        
	        for (j = 0; j < height; ++j)
	        {
	            y = (float)(j%gran)/gran;
	            yInv = y-1f;
	            gJ = j/gran;
	            gJ1 = gJ+1;
	            sy = s_curve(y);
	            
	            for (k = 0; k < depth; ++k)
	            {
	                z = (float)(k%gran)/gran;
	                zInv = z-1f;
	                gK = k/gran;
	                gK1 = gK+1;
	                sz = s_curve(z);
	                
	                gradient = gradients3D(gI,  gJ,  gK );
	                o = gradient[0]*x    + gradient[1]*y    + gradient[2]*z;
	                
	                gradient = gradients3D(gI1, gJ,  gK );
	                p = gradient[0]*xInv + gradient[1]*y    + gradient[2]*z;
	                
	                gradient = gradients3D(gI,  gJ1, gK );
	                q = gradient[0]*x    + gradient[1]*yInv + gradient[2]*z;
	                
	                gradient = gradients3D(gI1, gJ1, gK );
	                r = gradient[0]*xInv + gradient[1]*yInv + gradient[2]*z;
	                
	                gradient = gradients3D(gI,  gJ,  gK1);
	                s = gradient[0]*x    + gradient[1]*y    + gradient[2]*zInv;
	                
	                gradient = gradients3D(gI1, gJ,  gK1);
	                t = gradient[0]*xInv + gradient[1]*y    + gradient[2]*zInv;
	                
	                gradient = gradients3D(gI,  gJ1, gK1);
	                u = gradient[0]*x    + gradient[1]*yInv + gradient[2]*zInv;
	                
	                gradient = gradients3D(gI1, gJ1, gK1);
	                v = gradient[0]*xInv + gradient[1]*yInv + gradient[2]*zInv;
	                
	                
	                o =  o + sx*(p-o);
	                p =  q + sx*(r-q);
	                q =  s + sx*(t-s);
	                r =  u + sx*(v-u);
	                
	                o =  o + sy*(p-o);
	                p =  q + sy*(r-q);
	                
	                result[(int)n++] = o + sz*(p-o);
	            }
	        }
	    }
	    
	    return result;
	}*/
}
