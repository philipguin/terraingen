package src;

public final class Array2D
{
	public final int width, height;
	private final float[] values;
	
	public Array2D(float[] values, int width, int height)
	{
		this.width = width;
		this.height = height;
		this.values = new float[width * height];
	}
	
	public Array2D(int width, int height)
	{
		this(new float[width * height], width, height);
	}
	
	public final void set(int i, int j, float value)
	{
		values[j * width + i] = value;
	}
	
	public final float get(int i, int j)
	{
		return values[j * width + i];
	}
	
	public final void set1D(int index, float value)
	{
		values[index] = value;
	}

	public final float get1D(int index)
	{
		return values[index];
	}
	
	public final float[] exposeValues()
	{
		return values;
	}
}
