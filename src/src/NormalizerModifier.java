package src;

public class NormalizerModifier implements IModifier<Array2D>
{
	private final float min, range;
	
	public NormalizerModifier(float min, float max)
	{
		this.min = min;
		this.range = max - min;
	}
	
	public void modify(Array2D values)
	{
		if (values.width * values.height == 0)
			return;
		
	    float least = values.get1D(0), greatest = values.get1D(0);
	    
	    int end = values.width * values.height;
	    
	    for (int v = 1; v < end; ++v)
	    {
	    	if (values.get1D(v) < least)
	    		least = values.get1D(v);
	    	
	    	if (values.get1D(v) > greatest)
	    		greatest = values.get1D(v);
	    }
	    
	    final float oneOverCurrentRange = 1f / (greatest - least);
	    
	    for (int v = 0; v < values.width * values.height; ++v)
	    {
	    	values.set1D(v, (values.get1D(v) - least) * oneOverCurrentRange * range + min);
	    }
	    
	    //System.out.println("least: " + least + "\tgreatest: " + greatest);
	}
}
