package src;

public class Array3D
{
	public final int width, height, depth;
	private final float[] values;
	
	public Array3D(int width, int height, int depth)
	{
		this.width = width;
		this.height = height;
		this.depth = depth;
		this.values = new float[width * height * depth];
	}
	
	public final void set(int i, int j, int k, float value)
	{
		values[(k * height + j) * width + i] = value;
	}
	
	public final float get(int i, int j, int k)
	{
		return values[(k * height + j) * width + i];
	}
	
	public final float[] exposeValues()
	{
		return values;
	}
}
